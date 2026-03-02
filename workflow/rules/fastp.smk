rule fastp:
    input:
        unpack(get_fastqs)
    output:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        report="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
        signal=touch("results/fastp/{sample}.ok"),
    container: config["environments"]["custom"]
    log: "logs/fastp/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.report} -j {output.json} \
              -w {resources.threads}
        """
