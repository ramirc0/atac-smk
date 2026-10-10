# Library and signal QC metrics.


rule pbc_qc:
    input:
        f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
    output:
        f"{OUTDIR}/{{sample}}/pbc_qc/{{sample}}.pbc.json",
    log:
        f"{LOGDIR}/pbc_qc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/pbc_qc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/pbc_qc.py \
            --bam {input:q} \
            --output {output:q}
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


rule frip_all:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign.gz",
        peaks=f"{OUTDIR}/{{sample}}/macs3_callpeak/{{sample}}_peaks.narrowPeak",
    output:
        f"{OUTDIR}/{{sample}}/frip_all/{{sample}}.txt",
    log:
        f"{LOGDIR}/frip_all/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/frip_all/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/calculate_frip.py \
            --tagalign {input.tagalign:q} \
            --peaks {input.peaks:q} \
            --output {output:q}
        """
