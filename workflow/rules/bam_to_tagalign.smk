rule bam_to_tagalign:
    input: 
        bam = f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam"
    output: 
        bedpe = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.bedpe",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    container: config["container"]
    log: 
        f"{config['logdir']}/bedpe_tagalign/{{sample}}-{{genome}}.log"
    params:
        prefix = lambda wildcards: "-".join(wildcards),
    shell:
        """
        exec >> {log} 2>&1
        python ./workflow/rules/scripts/bam_to_tagalign.py {input.bam} {output.bedpe} {RESULTS_DIR}/bedpe_tagalign/{params.prefix}.tagalign --threads {resources.threads}
        """