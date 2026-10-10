rule bowtie2_align:
    input:
        index_signal=rules.bowtie2_build.output.signal,
        fastp_signal=rules.fastp.output.signal,
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
    output:
        bam="results/bowtie2_align/{sample}.bam",
        flagstat_json="results/bowtie2_align/{sample}.flagstat.json",
        signal=touch("results/bowtie2_align/{sample}.ok"),
    params:
        k=config['params']['k'],
        tmp_dir=config["input"]["tmp_dir"],
        index_prefix=rules.bowtie2_build.params.output_prefix,
    software: SOFTWARE_ENV
    log:
        "logs/bowtie2_align/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        bowtie2 -X 2000 --mm -k {params.k} \
            --threads {resources.threads} \
            -x {params.index_prefix} \
            -1 {input.r1} \
            -2 {input.r2} | \
        samtools view -@ {resources.threads} -S -1 -b /dev/stdin | \
        samtools sort -m 1536M -@ {resources.threads} -T {params.tmp_dir}/{wildcards.sample}.$$ -o {output.bam}
        samtools index -@ {resources.threads} {output.bam}
        samtools flagstat -@ {resources.threads} --output-fmt json {output.bam} > {output.flagstat_json}
        """
