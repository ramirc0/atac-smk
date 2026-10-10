rule macs3_pvalue_signal:
    input:
        tagalign_signal=rules.bam_to_tagalign.output.signal,
        macs3_signal=rules.macs3_signal.output.signal,
        bowtie2_build_signal=rules.bowtie2_build.output.signal,
        tagalign=rules.bam_to_tagalign.output.tagalign,
        treat_pileup=rules.macs3_signal.output.bedgraph,
        control_lambda=rules.macs3_signal.output.control_lambda,
        chromsizes=rules.bowtie2_build.output.chromsizes,
    output:
        pvalue_bedgraph=temp("results/macs3_pvalue_signal/{sample}_ppois.bdg"),
        pvalue_bedgraph_sorted="results/macs3_pvalue_signal/{sample}_ppois_sorted.bdg",
        pvalue_signal="results/macs3_pvalue_signal/{sample}_ppois.bw",
        signal=touch("results/macs3_pvalue_signal/{sample}.ok"),
    software: SOFTWARE_ENV
    log: "logs/macs3_pvalue_signal/{sample}.log"
    params:
        output_prefix="results/macs3_pvalue_signal/{sample}",
    shell:
        """
        exec >> {log} 2>&1

        SVAL=$(echo "scale=6; $(zcat {input.tagalign} | wc -l) / 1000000.0" | bc)

        macs3 bdgcmp -t {input.treat_pileup} \
            -c {input.control_lambda} \
            --o-prefix {params.output_prefix} -m ppois -S $SVAL

        bedtools slop -i {output.pvalue_bedgraph} -g {input.chromsizes} -b 0 | \
            bedClip stdin {input.chromsizes} /dev/stdout | \
            LC_COLLATE=C sort -k1,1 -k2,2n -S {resources.mem_mb}M | \
            awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' > {output.pvalue_bedgraph_sorted}

        bedGraphToBigWig {output.pvalue_bedgraph_sorted} {input.chromsizes} {output.pvalue_signal}
        """
