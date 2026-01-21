rule fastp:
    output:
        r1="results/fastp/{sample}_R1.fastq.gz",
        r2="results/fastp/{sample}_R2.fastq.gz",
        report="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
    params:
        sample=lambda wildcards: wildcards.sample,
        r1=lambda wildcards: READ_LOOKUP[wildcards.sample][0],
        r2=lambda wildcards: READ_LOOKUP[wildcards.sample][1],
    container: config["container"]
    log: f"{config['logdir']}/fastp/{{sample}}.log"
    shell:
        """
        exec >> {log} 2>&1
        echo "$(date --iso=minutes): Started fastp for sample: {params.sample}"

        fastp -i {params.r1} -I {params.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.report} -j {output.json} \
              -w {resources.threads}
        
        echo "$(date --iso=minutes): Completed fastp for sample: {params.sample}"
        """