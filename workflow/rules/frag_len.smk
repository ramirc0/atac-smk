rule frag_len:
    input:
        signal=rules.filter.output.signal,
        bam=rules.filter.output.bam,
    output: 
        png="results/frag_len/{sample}.png", 
        metrics="results/frag_len/{sample}.txt",
        signal=touch("results/frag_len/{sample}.ok"),
    software: SOFTWARE_ENV
    log: 
        "logs/frag_len/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/plot_fragment_length_distr.py {input.bam} {resources.threads} {output.png} {output.metrics}
        """

