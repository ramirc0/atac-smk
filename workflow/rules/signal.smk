# MACS3 signal tracks as bigWig.


rule macs3_signal:
    input:
        f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign.gz",
    output:
        bedgraph=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_treat_pileup.bdg",
        control_lambda=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_control_lambda.bdg",
    params:
        outdir=f"{OUTDIR}/{{sample}}/macs3_signal",
    log:
        f"{LOGDIR}/macs3_signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_signal/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 callpeak \
            -f BED \
            -t {input:q} \
            -n {wildcards.sample:q} \
            --outdir {params.outdir:q} \
            --nomodel \
            --shift -75 \
            --keep-dup all \
            --extsize 150 \
            -g hs \
            -q 0.01 \
            -B
        """


rule bedgraphtobigwig:
    input:
        bedgraph=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_treat_pileup.bdg",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
        chromsizes_bed=f"{REFDIR}/bowtie2_build/reference.chromsizes.bed",
    output:
        bedgraph=temp(f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}.bg"),
        bigwig=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}.bigWig",
    log:
        f"{LOGDIR}/bedgraphtobigwig/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bedgraphtobigwig/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect \
            -a {input.bedgraph:q} \
            -b {input.chromsizes_bed:q} \
            -sorted \
            > {output.bedgraph:q}
        bedGraphToBigWig {output.bedgraph:q} {input.chromsizes:q} {output.bigwig:q}
        """


rule macs3_fold_change_signal:
    input:
        treat_pileup=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_treat_pileup.bdg",
        control_lambda=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_control_lambda.bdg",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        bedgraph=temp(f"{OUTDIR}/{{sample}}/macs3_fold_change_signal/{{sample}}_FE.bdg"),
        bedgraph_sorted=f"{OUTDIR}/{{sample}}/macs3_fold_change_signal/{{sample}}_FE_sorted.bdg",
        bigwig=f"{OUTDIR}/{{sample}}/macs3_fold_change_signal/{{sample}}_FE.bw",
    params:
        prefix=f"{OUTDIR}/{{sample}}/macs3_fold_change_signal/{{sample}}",
    log:
        f"{LOGDIR}/macs3_fold_change_signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_fold_change_signal/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        macs3 bdgcmp \
            -t {input.treat_pileup:q} \
            -c {input.control_lambda:q} \
            --o-prefix {params.prefix:q} \
            -m FE

        bedtools slop -i {output.bedgraph:q} -g {input.chromsizes:q} -b 0 \
            | bedClip stdin {input.chromsizes:q} /dev/stdout \
            | LC_COLLATE=C sort -k1,1 -k2,2n -S {resources.mem_mb:q}M \
            | awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' \
            > {output.bedgraph_sorted:q}

        bedGraphToBigWig {output.bedgraph_sorted:q} {input.chromsizes:q} {output.bigwig:q}
        """


rule macs3_pvalue_signal:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign.gz",
        treat_pileup=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_treat_pileup.bdg",
        control_lambda=f"{OUTDIR}/{{sample}}/macs3_signal/{{sample}}_control_lambda.bdg",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        bedgraph=temp(f"{OUTDIR}/{{sample}}/macs3_pvalue_signal/{{sample}}_ppois.bdg"),
        bedgraph_sorted=f"{OUTDIR}/{{sample}}/macs3_pvalue_signal/{{sample}}_ppois_sorted.bdg",
        bigwig=f"{OUTDIR}/{{sample}}/macs3_pvalue_signal/{{sample}}_ppois.bw",
    params:
        prefix=f"{OUTDIR}/{{sample}}/macs3_pvalue_signal/{{sample}}",
    log:
        f"{LOGDIR}/macs3_pvalue_signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_pvalue_signal/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        SVAL=$(echo "scale=6; $(zcat {input.tagalign:q} | wc -l) / 1000000.0" | bc)

        macs3 bdgcmp \
            -t {input.treat_pileup:q} \
            -c {input.control_lambda:q} \
            --o-prefix {params.prefix:q} \
            -m ppois \
            -S $SVAL

        bedtools slop -i {output.bedgraph:q} -g {input.chromsizes:q} -b 0 \
            | bedClip stdin {input.chromsizes:q} /dev/stdout \
            | LC_COLLATE=C sort -k1,1 -k2,2n -S {resources.mem_mb:q}M \
            | awk 'BEGIN{{OFS="\t"}}{{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {{print $0}}; prev_chr=$1; prev_chr_e=$3;}}' \
            > {output.bedgraph_sorted:q}

        bedGraphToBigWig {output.bedgraph_sorted:q} {input.chromsizes:q} {output.bigwig:q}
        """
