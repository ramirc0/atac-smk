rule macs3_pseudoreps:
    input:
        signal=rules.pseudoreps.output.signal,
        tagalign="results/pseudoreps/{sample}_{pseudorep}.tagalign",
    output:
        narrowpeak="results/macs3_pseudoreps/{sample}_{pseudorep}_peaks.narrowPeak",
        signal=touch("results/macs3_pseudoreps/{sample}.{pseudorep}.ok"),
    params:
        parent_dir="results/macs3_pseudoreps/",
        p=config["params"]["p"],
    software: SOFTWARE_ENV
    log: "logs/macs3_pseudoreps/{sample}.{pseudorep}.log"
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {wildcards.sample}_{wildcards.pseudorep} --outdir {params.parent_dir} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs --pvalue {params.p} -B
        touch {output.signal}
        """
