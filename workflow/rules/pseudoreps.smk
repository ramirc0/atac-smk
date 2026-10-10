rule pseudoreps:
    input: 
        signal=rules.bam_to_tagalign.output.signal,
        bedpe=rules.bam_to_tagalign.output.bedpe,
    output: 
        tagalign=["results/pseudoreps/{sample}_pseudorep1.tagalign", "results/pseudoreps/{sample}_pseudorep2.tagalign"],
        signal=touch("results/pseudoreps/{sample}.ok"),
    software: SOFTWARE_ENV
    log: "logs/pseudoreps/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/bedpe_to_pseudoreps.py {input.bedpe} {output.tagalign[0]} {output.tagalign[1]}
        """

