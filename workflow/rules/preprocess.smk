# Reads to a filtered BAM and its tagAlign, as the ENCODE ATAC-seq pipeline does.


# Report only: every fastp trim and filter is off. Trimming is the trim rule's job.
rule fastp_qc:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        html=f"{OUTDIR}/{{sample}}/fastp/{{sample}}.html",
        json=f"{OUTDIR}/{{sample}}/fastp/{{sample}}.json",
    log:
        f"{LOGDIR}/fastp_qc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fastp_qc/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        fastp \
            -i {input.r1:q} \
            -I {input.r2:q} \
            --disable_adapter_trimming \
            --disable_quality_filtering \
            --disable_length_filtering \
            --disable_trim_poly_g \
            -h {output.html:q} \
            -j {output.json:q} \
            -w {threads}
        """


# ENCODE's detect_adapter.py picks the adapter seen most often in a mate's first
# 250,000 reads, ties in Illumina, Nextera, smallRNA order. Without an adapter in
# both mates, ENCODE leaves the reads untrimmed.
rule trim:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        r1=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R2.fastq.gz",
        json=f"{OUTDIR}/{{sample}}/trim/{{sample}}.cutadapt.json",
    log:
        f"{LOGDIR}/trim/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/trim/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        detect() {{
            awk 'NR % 4 == 2 {{
                    if (index($0, "AGATCGGAAGAGC")) illumina++
                    if (index($0, "CTGTCTCTTATA")) nextera++
                    if (index($0, "TGGAATTCTCGG")) smallrna++
                }}
                END {{
                    if (illumina > 0 && illumina >= nextera && illumina >= smallrna) print "AGATCGGAAGAGC"
                    else if (nextera > 0 && nextera >= smallrna) print "CTGTCTCTTATA"
                    else if (smallrna > 0) print "TGGAATTCTCGG"
                }}' <(zcat -f "$1" | head -n 1000000)
        }}
        adapter1=$(detect {input.r1:q})
        adapter2=$(detect {input.r2:q})
        echo "adapters: R1 '$adapter1', R2 '$adapter2'"

        if [[ -n "$adapter1" && -n "$adapter2" ]]; then
            cutadapt \
                -j {threads} \
                -e 0.1 \
                -m 5 \
                -a "$adapter1" \
                -A "$adapter2" \
                --json {output.json:q} \
                -o {output.r1:q} \
                -p {output.r2:q} \
                {input.r1:q} \
                {input.r2:q}
        else
            cp {input.r1:q} {output.r1:q}
            cp {input.r2:q} {output.r2:q}
            echo '{{}}' > {output.json:q}
        fi
        """


# Run-independent: built once under results/refs/, reused across runs.
# conda's bowtie2-build lacks SAIS, so its size check alone picks a small .bt2
# index for GRCh38. --large-index keeps the .bt2l index the pipeline expects.
rule bowtie2_build:
    input:
        REFERENCE_FASTA,
    output:
        index=BOWTIE2_INDEX_FILES,
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
        chromsizes_bed=f"{REFDIR}/bowtie2_build/reference.chromsizes.bed",
    params:
        prefix=BOWTIE2_INDEX,
    log:
        f"{LOGDIR}/bowtie2_build.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_build.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools faidx -@ {threads:q} {input:q} -o - | cut -f 1,2 > {output.chromsizes:q}
        cat {output.chromsizes:q} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,0,$2}}' \
            | sort -k1,1 -k2,2n \
            > {output.chromsizes_bed:q}
        bowtie2-build \
            --large-index \
            --threads {threads:q} \
            {input:q} \
            {params.prefix:q}
        """


# -k is one more than the alignments a read may have, so filter can tell which
# reads have too many.
rule bowtie2_align:
    input:
        index=BOWTIE2_INDEX_FILES,
        r1=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R2.fastq.gz",
    output:
        bam=f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
        flagstat=f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.flagstat.json",
    params:
        k=lambda w: config["params"]["multimapping"] + 1,
        index=BOWTIE2_INDEX,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/bowtie2_align/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_align/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        bowtie2 \
            -k {params.k:q} \
            -X 2000 \
            --mm \
            --threads {threads:q} \
            -x {params.index:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            | samtools sort \
                -m 1536M \
                -@ {threads:q} \
                -T "$tmp/sort" \
                -o {output.bam:q}
        samtools index -@ {threads:q} {output.bam:q}
        samtools flagstat -@ {threads:q} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


# ENCODE's paired-end filter with multimapping on. No MAPQ cutoff: reads with up to
# --multimapping alignments keep their primary one, and read names with more are
# dropped, as assign_multimappers.py does. That script also drops the last read
# name; this does not. -F 524 drops unmapped, mate unmapped and QC fail reads.
# -F 1804 adds secondary and duplicate. -f 2 keeps proper pairs.
rule filter:
    input:
        f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        bai=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam.bai",
        flagstat=f"{OUTDIR}/{{sample}}/filter/{{sample}}.flagstat.json",
        nmsrt=temp(f"{OUTDIR}/{{sample}}/filter/{{sample}}.nmsrt.bam"),
    params:
        max_lines=lambda w: 2 * config["params"]["multimapping"],
        mito=lambda w: config["mito_chr_name"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/filter/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/filter/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        samtools view -@ {threads:q} -u -F 524 -f 2 {input:q} \
            | samtools sort -@ {threads:q} -m 1536M -n -u -T "$tmp/nmsrt" - \
            | samtools view -h - \
            | awk -v max={params.max_lines:q} '
                /^@/ {{ print; next }}
                $1 != name {{ if (n <= max) printf "%s", group; name = $1; n = 0; group = "" }}
                {{ n++; group = group $0 "\n" }}
                END {{ if (n <= max) printf "%s", group }}' \
            | samtools fixmate -@ {threads:q} -r -u - - \
            | samtools view -@ {threads:q} -1 -F 1804 -f 2 -o {output.nmsrt:q} -
        samtools sort -@ {threads:q} -m 1536M -l 1 -T "$tmp/sort" -o "$tmp/filt.bam" {output.nmsrt:q}
        # ENCODE's sizes. They change which copy of a duplicate survives.
        sambamba markdup \
            -t {threads:q} \
            --tmpdir "$tmp" \
            --hash-table-size=17592186044416 \
            --overflow-list-size=20000000 \
            --io-buffer-size=256 \
            "$tmp/filt.bam" \
            "$tmp/dupmark.bam"
        samtools view \
            -@ {threads:q} \
            -b \
            -F 1804 \
            -f 2 \
            -e 'rname != "{params.mito}"' \
            -o {output.bam:q} \
            "$tmp/dupmark.bam"
        samtools index -@ {threads:q} {output.bam:q}
        samtools flagstat -@ {threads:q} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


# ENCODE's bam2ta: one line per mate, mate 1 first, in name order. Tn5 shift: +4 on
# the + strand, -5 on the - strand, clamped to at least 1 bp.
rule tagalign:
    input:
        f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
    output:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        bedpe=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.bedpe",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/tagalign/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/tagalign/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        samtools sort -@ {threads:q} -m 1536M -n -u -T "$tmp/nmsrt" {input:q} \
            | bedtools bamtobed -bedpe -mate1 -i stdin \
            | tee {output.bedpe:q} \
            | awk '
                function tag(chrom, start, end, strand) {{
                    if (strand == "+") start += 4
                    else if (strand == "-") end -= 5
                    if (start >= end) {{ if (strand == "+") start = end - 1; else end = start + 1 }}
                    print chrom "\t" start "\t" end "\tN\t1000\t" strand
                }}
                {{ tag($1, $2, $3, $9); tag($4, $5, $6, $10) }}' \
            | python -m zlib_ng.gzip_ng -cn \
            > {output.tagalign:q}
        """
