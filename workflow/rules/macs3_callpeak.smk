rule macs3_callpeak:
    input:
        signal=rules.bam_to_tagalign.output.signal,
        tagalign=rules.bam_to_tagalign.output.tagalign,
    output:
        narrowpeak="results/macs3_callpeak/{sample}_peaks.narrowPeak",
        signal=touch("results/macs3_callpeak/{sample}.ok"),
    software: SOFTWARE_ENV
    params:
        q=config["params"]["q"],
        parent_dir="results/macs3_callpeak/",
    log: 
        "logs/macs3_callpeak/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {wildcards.sample} --outdir {params.parent_dir} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q {params.q}
        """

