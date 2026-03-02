rule tss_enrichment:
    input: 
        signal=rules.filter.output.signal,
        bam=rules.filter.output.bam,
        tss_reference=config['resources']['tss_reference'],
    output:
        png="results/tss_enrichment/{sample}.png", 
        metrics="results/tss_enrichment/{sample}.txt",
        signal=touch("results/tss_enrichment/{sample}.ok"),
    container: config["environments"]["custom"]
    log: 
        "logs/tss_enrichment/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/tss_enrichment.py {input.bam} {input.tss_reference} {output.png} {output.metrics}
        """

