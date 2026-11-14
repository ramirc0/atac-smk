# The run's QC report: every sample's QC in one table, graded on ENCODE's ATAC
# portal audits (workflow/report/qc.yaml), with interactive figures.

import subprocess

REPORT_DIR = f"{OUTDIR}/report"
# The header names the pipeline commit, so a new commit renders the page again.
PIPELINE_VERSION = subprocess.run(
    ["git", "-C", workflow.basedir, "describe", "--always", "--dirty"],
    capture_output=True,
    text=True,
).stdout.strip()


rule gather_qc:
    input:
        expand(
            [
                f"{OUTDIR}/{{sample}}/fastp/{{sample}}.json",
                f"{OUTDIR}/{{sample}}/trim/{{sample}}.cutadapt.json",
                f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.flagstat.json",
                f"{OUTDIR}/{{sample}}/filter/{{sample}}.flagstat.json",
                f"{OUTDIR}/{{sample}}/pbc_qc/{{sample}}.pbc.json",
                f"{OUTDIR}/{{sample}}/frac_mito/{{sample}}.frac_mito.qc",
                f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.nucleosomal.qc",
                f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.tss_enrich.qc",
                f"{OUTDIR}/{{sample}}/frip/{{sample}}.tsv",
                f"{OUTDIR}/{{sample}}/annot_enrich/{{sample}}.annot_enrich.qc",
                f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.gc_summary.txt",
                f"{OUTDIR}/{{sample}}/jsd/{{sample}}.jsd.qc",
                f"{OUTDIR}/{{sample}}/peaks/{{sample}}.bfilt.narrowPeak.gz",
                f"{OUTDIR}/{{sample}}/overlap/{{sample}}.overlap.bfilt.narrowPeak.gz",
                f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.bfilt.narrowPeak.gz",
            ],
            sample=SAMPLES,
        ),
    output:
        f"{REPORT_DIR}/qc_metrics.tsv",
    params:
        run=OUTDIR,
        samples=SAMPLES,
    log:
        f"{LOGDIR}/gather_qc/all.txt",
    benchmark:
        f"{BENCHDIR}/gather_qc/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/gather_qc.py \
            --run-dir {params.run:q} \
            --samples {params.samples:q} \
            --threads {threads:q} \
            --output {output:q}

        sync {output:q}
        """


rule qc_figures:
    input:
        expand(
            [
                f"{OUTDIR}/{{sample}}/fastp/{{sample}}.json",
                f"{OUTDIR}/{{sample}}/tss_enrichment/{{sample}}.tss_profile.tsv",
                f"{OUTDIR}/{{sample}}/frag_len/{{sample}}.insert_size.txt",
                f"{OUTDIR}/{{sample}}/gc_bias/{{sample}}.gc.txt",
                f"{OUTDIR}/{{sample}}/jsd/{{sample}}.fingerprint.tsv",
                f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.txt.gz",
            ],
            sample=SAMPLES,
        ),
    output:
        expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    params:
        run=OUTDIR,
        report_dir=REPORT_DIR,
        samples=SAMPLES,
        idr_threshold=lambda w: config["params"]["idr_thresh"],
    log:
        f"{LOGDIR}/qc_figures/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_figures/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_figures.py \
            --run-dir {params.run:q} \
            --report-dir {params.report_dir:q} \
            --samples {params.samples:q} \
            --idr-threshold {params.idr_threshold:q} \
            --threads {threads:q}

        sync {output:q}
        """


rule qc_report:
    input:
        spec=f"{workflow.basedir}/report/qc.yaml",
        metrics=f"{REPORT_DIR}/qc_metrics.tsv",
        figures=expand(f"{REPORT_DIR}/samples/{{sample}}.js", sample=SAMPLES),
    output:
        f"{REPORT_DIR}/qc_report.html",
    params:
        run=OUTDIR,
        run_id=RUN_ID,
        version=f"atac-smk {PIPELINE_VERSION}",
    log:
        f"{LOGDIR}/qc_report/all.txt",
    benchmark:
        f"{BENCHDIR}/qc_report/all.tsv"
    software:
        REPORT_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        qc-report \
            --spec {input.spec:q} \
            --metrics {input.metrics:q} \
            --run-dir {params.run:q} \
            --run-id {params.run_id:q} \
            --version {params.version:q} \
            --output {output:q}

        sync {output:q}
        """
