rule macs3_fold_change_signal:
    input:
        macs3_signal=rules.macs3_signal.output.signal,
        bowtie2_build_signal=rules.bowtie2_build.output.signal,
        treat_pileup=rules.macs3_signal.output.bedgraph,
        control_lambda=rules.macs3_signal.output.control_lambda,
        chromsizes=rules.bowtie2_build.output.chromsizes,
    output:
        fc_bedgraph=temp("results/macs3_fold_change_signal/{sample}_FE.bdg"),
        fc_bedgraph_sorted="results/macs3_fold_change_signal/{sample}_FE_sorted.bdg",
        fc_signal="results/macs3_fold_change_signal/{sample}_FE.bw",
        signal=touch("results/macs3_fold_change_signal/{sample}.ok"),
    software: SOFTWARE_ENV
    log: "logs/macs3_fold_change_signal/{sample}.log"
    params:
        output_prefix="results/macs3_fold_change_signal/{sample}",
    shell:
        """
        exec >> {log} 2>&1

        macs3 bdgcmp -t {input.treat_pileup} \
            -c {input.control_lambda} \
            --o-prefix {params.output_prefix} -m FE

        bedtools slop -i {output.fc_bedgraph} -g {input.chromsizes} -b 0 | \
            bedClip stdin {input.chromsizes} /dev/stdout | \
            LC_COLLATE=C sort -k1,1 -k2,2n -S {resources.mem_mb}M | \
            awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' > {output.fc_bedgraph_sorted}

        bedGraphToBigWig {output.fc_bedgraph_sorted} {input.chromsizes} {output.fc_signal}
        """
