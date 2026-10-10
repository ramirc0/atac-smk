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
        """
