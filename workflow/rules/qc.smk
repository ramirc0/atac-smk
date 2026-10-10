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


rule frag_len:
    input:
        f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
    output:
        plot=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.png",
        metrics=f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.txt",
    log:
        f"{LOGDIR}/frag_len/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frag_len/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_fragment_length_distr.py \
            --bam {input:q} \
            --plot {output.plot:q} \
            --output {output.metrics:q}
        """


rule tss_enrichment:
    input:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        tss=TSS_REFERENCE,
    output:
        plot=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.svg",
        plot_png=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.png",
        metrics=f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.txt",
    log:
        f"{LOGDIR}/tss_enrichment/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/tss_enrichment/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/tss_enrichment.py \
            --bam {input.bam:q} \
            --tss {input.tss:q} \
            --plot {output.plot:q} \
            --output {output.metrics:q}
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
