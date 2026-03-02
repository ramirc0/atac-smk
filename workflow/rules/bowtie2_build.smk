rule bowtie2_build:
    input:
        fa=config["resources"]["reference_fasta"]
    output:
        index = multiext(
            "results/bowtie2_build/reference",
            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l",
            ".rev.1.bt2l", ".rev.2.bt2l"
        ),
        chromsizes="results/bowtie2_build/reference.chromsizes",
        chromsizes_bed="results/bowtie2_build/reference.chromsizes.bed",
        signal=touch("results/bowtie2_build/reference.ok"),
    params:
        output_prefix = "results/bowtie2_build/reference",
    container: config["environments"]["custom"]
    log: "logs/bowtie2_build.log"
    shell:
        """
        exec >> {log} 2>&1
        samtools faidx -@ {resources.threads} {input.fa} -o - | cut -f 1,2 > {output.chromsizes}
        cat {output.chromsizes} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,0,$2}}' | sort -k1,1 -k2,2n > {output.chromsizes_bed}
        bowtie2-build --threads {resources.threads} {input.fa} {params.output_prefix}
        """

