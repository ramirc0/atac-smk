# Library and signal QC metrics.


# ENCODE's library complexity, on the filtered reads before duplicate removal.
rule pbc_qc:
    input:
        f"{OUTDIR}/{{sample}}/filter/{{sample}}.nmsrt.bam",
    output:
        f"{OUTDIR}/{{sample}}/pbc_qc/{{sample}}.pbc.json",
    params:
        mito=lambda w: config["mito_chr_name"],
    log:
        f"{LOGDIR}/pbc_qc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/pbc_qc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools bamtobed -bedpe -i {input:q} \
            | python workflow/scripts/pbc_qc.py --mito {params.mito:q} --output {output:q}

        sync {output:q}
        """


# ENCODE's fragment length QC, on every record. ENCODE stops Picard after the
# first 5M records, which in a sorted BAM come mostly from chr1.
rule frag_len:
    input:
        f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
    output:
        qc=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.nucleosomal.qc",
        plot=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.png",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/frag_len/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frag_len/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        picard -Xmx4g CollectInsertSizeMetrics \
            -I {input:q} \
            -O "$tmp/hist.txt" \
            -H "$tmp/hist.pdf" \
            -W 1000 \
            -VERBOSITY ERROR \
            -QUIET true
        python workflow/scripts/fraglen_qc.py \
            --histogram "$tmp/hist.txt" \
            --qc {output.qc:q} \
            --plot {output.plot:q}

        sync {output:q}
        """


rule tss_enrichment:
    input:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        bai=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam.bai",
        tss=TSS_REFERENCE,
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
        r1=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R1.fastq.gz",
    output:
        qc=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.tss_enrich.qc",
        plot=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.png",
    log:
        f"{LOGDIR}/tss_enrichment/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/tss_enrichment/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/tss_enrichment.py \
            --bam {input.bam:q} \
            --tss {input.tss:q} \
            --chromsizes {input.chromsizes:q} \
            --fastq {input.r1:q} \
            --qc {output.qc:q} \
            --plot {output.plot:q}

        sync {output:q}
        """


# ENCODE's FRiP: tagAlign lines that touch a peak, for each filtered peak set.
rule frip:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        peaks=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.bfilt.narrowPeak.gz",
        idr=f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.bfilt.narrowPeak.gz",
        overlap=f"{OUTDIR}/{{sample}}/overlap/{{sample}}.overlap.bfilt.narrowPeak.gz",
    output:
        f"{OUTDIR}/{{sample}}/frip/{{sample}}.tsv",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/frip/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frip/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        # ENCODE decompresses first; bedtools mishandles two gzipped inputs.
        python -m zlib_ng.gzip_ng -cd {input.tagalign:q} > "$tmp/ta"
        reads=$(wc -l < "$tmp/ta")
        frip() {{
            python -m zlib_ng.gzip_ng -cd "$2" > "$tmp/peaks"
            hits=$(bedtools intersect -nonamecheck -a "$tmp/ta" -b "$tmp/peaks" -wa -u | wc -l)
            awk -v set="$1" -v hits="$hits" -v reads="$reads" \
                'BEGIN {{ printf "%s\t%d\t%d\t%.6f\n", set, hits, reads, hits / reads }}'
        }}
        {{
            printf "peak_set\treads_in_peaks\treads\tfrip\n"
            frip peaks {input.peaks:q}
            frip idr {input.idr:q}
            frip overlap {input.overlap:q}
        }} > {output:q}

        sync {output:q}
        """


# ENCODE's frac_mito: reads that align to the mitochondrial contig alone, over those
# plus the reads mapped elsewhere in the main alignment. A read counts once, as
# mapped if any of its alignments is, as SAMstats counts it.
rule frac_mito:
    input:
        index=MITO_INDEX_FILES,
        r1=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/trim/{{sample}}_R2.fastq.gz",
        bam=f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
    output:
        f"{OUTDIR}/{{sample}}/frac_mito/{{sample}}.frac_mito.qc",
    params:
        k=lambda w: config["params"]["multimapping"] + 1,
        index=MITO_INDEX,
        mito=lambda w: config["mito_chr_name"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/frac_mito/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frac_mito/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        mito=$(bowtie2 \
            -k {params.k:q} \
            -X 2000 \
            --mm \
            --threads {threads:q} \
            -x {params.index:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            | samtools view -c -F 0x904 -)
        # Read name plus mate (64 = first, 128 = second), from every mapped alignment.
        non_mito=$(samtools view -@ {threads:q} -F 0x804 -e 'rname != "{params.mito}"' {input.bam:q} \
            | awk '{{ print $1 "\t" int($2 / 64) % 4 }}' \
            | sort -u -S 2G --parallel {threads:q} -T "$tmp" \
            | wc -l)
        awk -v n="$non_mito" -v m="$mito" 'BEGIN {{
            printf "non_mito_reads\t%d\nmito_reads\t%d\nfrac_mito_reads\t%.6f\n", n, m, m / (n + m)
        }}' > {output:q}

        sync {output:q}
        """


# ENCODE's annot_enrich: tagAlign lines in each merged region set.
rule annot_enrich:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        dhs=DHS,
        blacklist=BLACKLIST,
        promoters=PROMOTERS,
        enhancers=ENHANCERS,
    output:
        f"{OUTDIR}/{{sample}}/annot_enrich/{{sample}}.annot_enrich.qc",
    log:
        f"{LOGDIR}/annot_enrich/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/annot_enrich/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        reads=$(python -m zlib_ng.gzip_ng -cd {input.tagalign:q} | wc -l)
        fraction() {{
            hits=$(bedtools sort -i "$2" \
                | bedtools merge -i stdin \
                | bedtools intersect -u -nonamecheck -a {input.tagalign:q} -b stdin \
                | wc -l)
            awk -v name="$1" -v hits="$hits" -v reads="$reads" \
                'BEGIN {{ printf "%s\t%d\t%.6f\n", name, hits, hits / reads }}'
        }}
        {{
            fraction fraction_of_reads_in_universal_DHS_regions {input.dhs:q}
            fraction fraction_of_reads_in_blacklist_regions {input.blacklist:q}
            fraction fraction_of_reads_in_promoter_regions {input.promoters:q}
            fraction fraction_of_reads_in_enhancer_regions {input.enhancers:q}
        }} > {output:q}

        sync {output:q}
        """


# ENCODE's GC bias: Picard's metrics, replotted.
rule gc_bias:
    input:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        fasta=REFERENCE_FASTA,
    output:
        metrics=f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.gc.txt",
        summary=f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.gc_summary.txt",
        plot=f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.png",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/gc_bias/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/gc_bias/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        picard -Xmx4g CollectGcBiasMetrics \
            -R {input.fasta:q} \
            -I {input.bam:q} \
            -O {output.metrics:q} \
            -CHART "$tmp/gc.pdf" \
            -S {output.summary:q} \
            -ASSUME_SORTED false \
            -VERBOSITY ERROR \
            -QUIET true
        python workflow/scripts/plot_gc.py --metrics {output.metrics:q} --plot {output.plot:q}

        sync {output:q}
        """


# ENCODE's JSD: deeptools fingerprint metrics on the blacklist-filtered BAM.
rule jsd:
    input:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        blacklist=BLACKLIST,
    output:
        qc=f"{OUTDIR}/{{sample}}/jsd/{{sample}}.jsd.qc",
        plot=f"{OUTDIR}/{{sample}}/jsd/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/jsd/{{sample}}.png",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/jsd/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/jsd/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        bedtools intersect -nonamecheck -v -abam {input.bam:q} -b {input.blacklist:q} > "$tmp/bfilt.bam"
        samtools index "$tmp/bfilt.bam"
        plotFingerprint \
            -b "$tmp/bfilt.bam" \
            --labels {wildcards.sample:q} \
            --outQualityMetrics {output.qc:q} \
            --outRawCounts "$tmp/counts.tab" \
            --minMappingQuality 30 \
            --numberOfProcessors {threads:q} \
            --plotFile "$tmp/fingerprint.png"
        python workflow/scripts/plot_fingerprint.py \
            --counts "$tmp/counts.tab" \
            --title {wildcards.sample:q} \
            --plot {output.plot:q}

        sync {output:q}
        """
