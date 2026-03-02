rule macs3_signal:
    input:
        signal=rules.bam_to_tagalign.output.signal,
        tagalign=rules.bam_to_tagalign.output.tagalign,
    output:
        bedgraph="results/macs3_signal/{sample}_treat_pileup.bdg",
        control_lambda="results/macs3_signal/{sample}_control_lambda.bdg",
        signal=touch("results/macs3_signal/{sample}.ok"),
    params:
        parent_dir="results/macs3_signal/",
    container: config["environments"]["custom"]
    log: "logs/macs3_signal/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {wildcards.sample} --outdir {params.parent_dir} --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q 0.01 -B
        """
