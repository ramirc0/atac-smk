rule bam_to_tagalign:
    input: 
        signal=rules.filter.output.signal,
        bam=rules.filter.output.bam,
    output: 
        bedpe="results/bedpe_tagalign/{sample}.bedpe",
        tagalign="results/bedpe_tagalign/{sample}.tagalign.gz",
        signal=touch("results/bedpe_tagalign/{sample}.ok"),
    params:
      tagalign_prefix="results/bedpe_tagalign/{sample}.tagalign",
    software: SOFTWARE_ENV
    log: "logs/bedpe_tagalign/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/bam_to_tagalign.py {input.bam} {output.bedpe} {params.tagalign_prefix} --threads {resources.threads}
        """
