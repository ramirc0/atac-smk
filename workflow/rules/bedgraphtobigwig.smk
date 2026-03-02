rule bedgraphtobigwig:
    input: 
        macs3_signal=rules.macs3_signal.output.signal,
        bowtie2_build_signal=rules.bowtie2_build.output.signal,
        bedgraph=rules.macs3_signal.output.bedgraph,
        chromsizes=rules.bowtie2_build.output.chromsizes,
        chromsizes_bed=rules.bowtie2_build.output.chromsizes_bed,
    output:
        bedgraph=temp("results/macs3_signal/{sample}.bg"),
        bigwig="results/macs3_signal/{sample}.bigWig",
        signal=touch("results/macs3_signal/{sample}.bedgraphtobigwig.ok"),
    container: config["environments"]["custom"]
    log: "logs/bedgraphtobigwig/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        bedtools intersect -a {input.bedgraph} -b {input.chromsizes_bed} -sorted > {output.bedgraph}
        bedGraphToBigWig {output.bedgraph} {input.chromsizes} {output.bigwig}
        """
