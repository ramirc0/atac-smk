# | Bit | Value | Excludes |
# |-----|-------|----------|
# | 4 | 4 | Unmapped |
# | 8 | 8 | Mate unmapped |
# | 256 | 256 | Secondary alignment |
# | 512 | 512 | Failed QC |
# | 1024 | 1024 | PCR duplicate |

rule filter:
    input: 
        bam="results/bowtie2_align/{sample}.bam",
        signal=rules.bowtie2_align.output.signal,
    output: 
        bam="results/filter/{sample}.bam",
        flagstat_json="results/filter/{sample}.flagstat.json",
        signal=touch("results/filter/{sample}.ok"),
    log: 
        "results/filter/{sample}.log"
    params:
        exclude_flag=config['params']['samtools_exclude_flag'],
        include_flag=config['params']['samtools_include_flag'],
    software: SOFTWARE_ENV
    shell:
        """
        exec >> {log} 2>&1

        sambamba markdup -t {resources.threads} {input.bam} /dev/stdout | \
        samtools view -@ {resources.threads} -F {params.exclude_flag} -f {params.include_flag} -q 30 -h -u /dev/stdin | \
        samtools sort -m 1536M -@ {resources.threads} -n /dev/stdin | \
        mbuffer -m 1G | \
        samtools fixmate -@ {resources.threads} -r /dev/stdin /dev/stdout | \
        samtools view -@ {resources.threads} -F {params.exclude_flag} -f {params.include_flag} -h -u /dev/stdin | \
        samtools sort -m 1536M -@ {resources.threads} /dev/stdin | \
        mbuffer -m 1G | \
        samtools view -@ {resources.threads} -h | \
        grep -v "chrM" | \
        samtools view -@ {resources.threads} -b > {output.bam}

        samtools quickcheck -v {output.bam}
        samtools index -@ {resources.threads} {output.bam}
        samtools flagstat -@ {resources.threads} --output-fmt json {output.bam} > {output.flagstat_json}
        """
