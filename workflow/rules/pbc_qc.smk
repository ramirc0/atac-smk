rule pbc_qc:
    input:
        signal=rules.bowtie2_align.output.signal,
        bam=rules.bowtie2_align.output.bam,
    output:
        pbc_json="results/pbc_qc/{sample}.pbc.json",
        signal=touch("results/pbc_qc/{sample}.ok"),
    software: SOFTWARE_ENV
    log: "logs/pbc_qc/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/scripts/pbc_qc.py --bam {input.bam} --output {output.pbc_json}
        """
