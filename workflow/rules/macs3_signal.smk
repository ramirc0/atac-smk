import math

rule macs3_signal:
    input:
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz",
    output:
        bedgraph = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg",
        control_lambda = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_control_lambda.bdg",
    container: config["container"]
    log: 
        f"{config['logdir']}/macs3_signal/{{sample}}-{{genome}}.log"
    params:
        prefix = lambda wildcards: "-".join(wildcards),
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {params.prefix} --outdir {RESULTS_DIR}/macs3_signal/ --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q .01 -B
        """

rule macs3_fold_change_signal:
    input:
        treat_pileup = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg",
        control_lambda = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_control_lambda.bdg",
        chromsizes = rules.bowtie2_build.output.chromsizes,
    output:
        fc_bedgraph = temp(f"{RESULTS_DIR}/macs3_fold_change_signal/{{sample}}-{{genome}}_FE.bdg"),
        fc_bedgraph_sorted = f"{RESULTS_DIR}/macs3_fold_change_signal/{{sample}}-{{genome}}_FE_sorted.bdg",
        fc_signal = f"{RESULTS_DIR}/macs3_fold_change_signal/{{sample}}-{{genome}}_FE.bw"
    container: config["container"]
    log: f"{config['logdir']}/macs3_fold_change_signal/{{sample}}-{{genome}}.log"
    params:
        output_prefix = f"{RESULTS_DIR}/macs3_fold_change_signal/{{sample}}-{{genome}}",
    shell:
        """
        exec >> {log} 2>&1

        MEM_MB=$((({resources.mem_mb} * 50) / 100))
        macs3 bdgcmp -t {input.treat_pileup} \
            -c {input.control_lambda} \
            --o-prefix {params.output_prefix} -m FE

        bedtools slop -i {output.fc_bedgraph} -g {input.chromsizes} -b 0 | \
            bedClip stdin {input.chromsizes} /dev/stdout | \
            LC_COLLATE=C sort -k1,1 -k2,2n -S ${{MEM_MB}}M | \
            awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' > {output.fc_bedgraph_sorted}

        bedGraphToBigWig {output.fc_bedgraph_sorted} {input.chromsizes} {output.fc_signal}
        """

rule macs3_pvalue_signal:
    input:
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz",
        treat_pileup = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_treat_pileup.bdg",
        control_lambda = f"{RESULTS_DIR}/macs3_signal/{{sample}}-{{genome}}_control_lambda.bdg",
        chromsizes = rules.bowtie2_build.output.chromsizes,
    output:
        pvalue_bedgraph = temp(f"{RESULTS_DIR}/macs3_pvalue_signal/{{sample}}-{{genome}}_ppois.bdg"),
        pvalue_bedgraph_sorted = f"{RESULTS_DIR}/macs3_pvalue_signal/{{sample}}-{{genome}}_ppois_sorted.bdg",
        pvalue_signal = f"{RESULTS_DIR}/macs3_pvalue_signal/{{sample}}-{{genome}}_ppois.bw"
    container: config["container"]
    log: f"{config['logdir']}/macs3_pvalue_signal/{{sample}}-{{genome}}.log"
    params:
        output_prefix = f"{RESULTS_DIR}/macs3_pvalue_signal/{{sample}}-{{genome}}",
    shell:
        """
        exec >> {log} 2>&1

        MEM_MB=$((({resources.mem_mb} * 50) / 100))
        SVAL=$(echo "scale=6; $(zcat {input.tagalign} | wc -l) / 1000000.0" | bc)

        macs3 bdgcmp -t {input.treat_pileup} \
            -c {input.control_lambda} \
            --o-prefix {params.output_prefix} -m ppois -S $SVAL

        bedtools slop -i {output.pvalue_bedgraph} -g {input.chromsizes} -b 0 | \
            bedClip stdin {input.chromsizes} /dev/stdout | \
            LC_COLLATE=C sort -k1,1 -k2,2n -S ${{MEM_MB}}M | \
            awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' > {output.pvalue_bedgraph_sorted}

        bedGraphToBigWig {output.pvalue_bedgraph_sorted} {input.chromsizes} {output.pvalue_signal}
        """

