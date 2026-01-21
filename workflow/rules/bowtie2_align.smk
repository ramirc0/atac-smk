rule bowtie2_align:
    input:
        index_dir = rules.bowtie2_build.output.index_dir,
        index = multiext(
            f"{RESULTS_DIR}/bowtie2_build/{{genome}}/{{genome}}.fa",
            ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l",
            ".rev.1.bt2l", ".rev.2.bt2l"
        ),
        R1=f"{RESULTS_DIR}/fastp/{{sample}}_R1.fastq.gz",
        R2=f"{RESULTS_DIR}/fastp/{{sample}}_R2.fastq.gz"
    output:
        bam = f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam",
        flagstat_json = f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.flagstat.json"
    params:
        k=config['k']
    container: config["container"]
    log:
        f"{config['logdir']}/bowtie2_align/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "Running bowtie2 align"
        bowtie2 -X 2000 --mm -k {params.k} \
            --threads {resources.threads} \
            -x {input.index_dir}/{wildcards.genome}.fa \
            --rg-id {wildcards.sample}-{wildcards.genome} \
            --rg SM:{wildcards.sample}-{wildcards.genome} \
            -1 {input.R1} \
            -2 {input.R2} | \
        samtools view -@ {resources.threads} -S -1 -b /dev/stdin | \
        samtools sort -@ {resources.threads} -T {wildcards.sample}-{wildcards.genome} -o {output.bam}
        samtools index -@ {resources.threads} {output.bam}
        samtools flagstat -@ {resources.threads} --output-fmt json {output.bam} > {output.flagstat_json}
        """

rule pbc_qc:
    input:
        bam = f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam",
    output:
        pbc_json = f"{RESULTS_DIR}/pbc_qc/{{sample}}-{{genome}}.pbc.json",
    container: config["container"]
    log: f"{config['logdir']}/pbc_qc/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python ./workflow/rules/scripts/pbc_qc.py {input.bam} {output.pbc_json}
        """
