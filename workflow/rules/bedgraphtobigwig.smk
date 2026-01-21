rule bedgraphtobigwig:
    input: 
        bedgraph = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg",
        chromsizes = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.chromsizes",
        chromsizes_bed = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.chromsizes.bed",
    output:
        bigwig = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}.bigWig"
    container: config["container"]
    log:
        f"{config['logdir']}/bedgraphtobigwig/{{sample}}-{{genome}}.log"
    params:
        prefix = lambda wildcards: "-".join(wildcards),
    shell:
        """
        exec >> {log} 2>&1
        echo "Fixing bedGraph"
        bedtools intersect -a {input.bedgraph} -b {input.chromsizes_bed} -sorted > {RESULTS_DIR}/macs3_signal/{params.prefix}.bg

        echo "Running bigtools bedgraphtobigwig"
        bedGraphToBigWig {RESULTS_DIR}/macs3_signal/{params.prefix}.bg {input.chromsizes} {output.bigwig}
        """