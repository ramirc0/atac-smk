def get_fasta(wildcards):
    genome = wildcards.genome 
    result = config['genomes'].get(genome)
    return result

rule bowtie2_build:
    input:
        fa = get_fasta
    output:
        index = multiext(
            f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.fa",
            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l",
            ".rev.1.bt2l", ".rev.2.bt2l"
        ),
        fa = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.fa",
        chromsizes = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.chromsizes",
        chromsizes_bed = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.chromsizes.bed",
        index_dir = directory(f"{RESULTS_DIR}/bowtie2_build/{{genome}}"),
        signal = f"{RESULTS_DIR}/bowtie2_build/{{genome}}/.continue"
    container: config["container"]
    log: f"{config['logdir']}/bowtie2_build/{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "Building Bowtie2 index for {wildcards.genome}"

        cp {input.fa} {output.fa}
        samtools faidx -@ {resources.threads} {input.fa} -o - | cut -f 1,2 > {output.chromsizes}
        cat {output.chromsizes} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,0,$2}}' | sort -k1,1 -k2,2n > {output.chromsizes_bed}

        bowtie2-build --threads {resources.threads} {input.fa} {RESULTS_DIR}/bowtie2_build/{wildcards.genome}/{wildcards.genome}.fa
        echo "Bowtie2 index built for {wildcards.genome}"
        touch {output.signal}
        """

