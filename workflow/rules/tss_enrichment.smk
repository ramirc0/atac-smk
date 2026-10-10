rule tss_enrichment:
    input: 
        signal=rules.filter.output.signal,
        bam=rules.filter.output.bam,
        tss_reference=config['resources']['tss_reference'],
    output:
        plot="results/tss_enrichment/{sample}.svg",
        plot_png="results/tss_enrichment/{sample}.png",
        metrics="results/tss_enrichment/{sample}.txt",
        signal=touch("results/tss_enrichment/{sample}.ok"),
    software: SOFTWARE_ENV
    log: 
        "logs/tss_enrichment/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/scripts/tss_enrichment.py --bam {input.bam} --tss {input.tss_reference} --plot {output.plot} --output {output.metrics}
        """

