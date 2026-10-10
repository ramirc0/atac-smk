# Peaks for the full sample and its pseudoreplicates, and their reproducibility.


rule macs3_callpeak:
    input:
        f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign.gz",
    output:
        f"{OUTDIR}/{{sample}}/macs3_callpeak/{{sample}}_peaks.narrowPeak",
    params:
        q=lambda w: config["params"]["q"],
        outdir=f"{OUTDIR}/{{sample}}/macs3_callpeak",
    log:
        f"{LOGDIR}/macs3_callpeak/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_callpeak/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 callpeak \
            -f BED \
            -t {input:q} \
            -n {wildcards.sample:q} \
            --outdir {params.outdir:q} \
            --nomodel \
            --shift -75 \
            --keep-dup all \
            --extsize 150 \
            -g hs \
            -q {params.q:q}
        """


rule pseudoreps:
    input:
        f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.bedpe",
    output:
        pseudorep1=f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}_pseudorep1.tagalign",
        pseudorep2=f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}_pseudorep2.tagalign",
    log:
        f"{LOGDIR}/pseudoreps/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/pseudoreps/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/bedpe_to_pseudoreps.py \
            --bedpe {input:q} \
            --pseudorep1 {output.pseudorep1:q} \
            --pseudorep2 {output.pseudorep2:q}
        """


rule macs3_pseudoreps:
    input:
        f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}_{{pseudorep}}.tagalign",
    output:
        f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_{{pseudorep}}_peaks.narrowPeak",
    params:
        p=lambda w: config["params"]["p"],
        outdir=f"{OUTDIR}/{{sample}}/macs3_pseudoreps",
    log:
        f"{LOGDIR}/macs3_pseudoreps/{{sample}}.{{pseudorep}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_pseudoreps/{{sample}}.{{pseudorep}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 callpeak \
            -f BED \
            -t {input:q} \
            -n {wildcards.sample:q}_{wildcards.pseudorep:q} \
            --outdir {params.outdir:q} \
            --nomodel \
            --shift -75 \
            --keep-dup all \
            --extsize 150 \
            -g hs \
            --pvalue {params.p:q} \
            -B
        """


rule idr:
    input:
        pseudorep1=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep1_peaks.narrowPeak",
        pseudorep2=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep2_peaks.narrowPeak",
    output:
        bed=f"{OUTDIR}/{{sample}}/idr/{{sample}}.IDR.bed",
        png=f"{OUTDIR}/{{sample}}/idr/{{sample}}.IDR.png",
    params:
        prefix=lambda w: f"{TMP_DIR}/{w.sample}",
    log:
        f"{LOGDIR}/idr/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/idr/{{sample}}.tsv",
    software:
        IDR_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        sort -k8,8nr {input.pseudorep1:q} > {params.prefix:q}_pseudorep1.sorted
        sort -k8,8nr {input.pseudorep2:q} > {params.prefix:q}_pseudorep2.sorted
        idr \
            --samples {params.prefix:q}_pseudorep1.sorted {params.prefix:q}_pseudorep2.sorted \
            --input-file-type narrowPeak \
            --rank p.value \
            --output-file {params.prefix:q}.IDR.unsorted.bed \
            --plot
        cat {params.prefix:q}.IDR.unsorted.bed | sort -k1,1 -k2,2n > {output.bed:q}
        mv {params.prefix:q}.IDR.unsorted.bed.png {output.png:q}
        rm {params.prefix:q}*
        """


rule overlap_peaks:
    input:
        pseudorep1=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep1_peaks.narrowPeak",
        pseudorep2=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep2_peaks.narrowPeak",
    output:
        f"{OUTDIR}/{{sample}}/overlap_peaks/{{sample}}_overlap_peaks.bed",
    log:
        f"{LOGDIR}/overlap_peaks/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/overlap_peaks/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect -a {input.pseudorep1:q} -b {input.pseudorep2:q} -u \
            | sort -k1,1 -k2,2n \
            > {output:q}
        """


rule overlap_peaks_abc:
    input:
        peaks=f"{OUTDIR}/{{sample}}/macs3_callpeak/{{sample}}_peaks.narrowPeak",
        pseudorep1=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep1_peaks.narrowPeak",
        pseudorep2=f"{OUTDIR}/{{sample}}/macs3_pseudoreps/{{sample}}_pseudorep2_peaks.narrowPeak",
    output:
        f"{OUTDIR}/{{sample}}/overlap_peaks_abc/{{sample}}.overlap_peaks.bed",
    log:
        f"{LOGDIR}/overlap_peaks_abc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/overlap_peaks_abc/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/overlap_peaks_abc.py \
            --peaks {input.peaks:q} \
            --pseudorep1 {input.pseudorep1:q} \
            --pseudorep2 {input.pseudorep2:q} \
            --output {output:q}
        """
