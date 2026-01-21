# def get_exclude_flag(wildcards):
#     if int(config['k']) == -1 or int(config['k']) > 1:
#         return("1548")
#     else:
#         return("1804")

# | Bit | Value | Excludes |
# |-----|-------|----------|
# | 4 | 4 | Unmapped |
# | 8 | 8 | Mate unmapped |
# | 256 | 256 | Secondary alignment |
# | 512 | 512 | Failed QC |
# | 1024 | 1024 | PCR duplicate |

rule filter:
    input: 
        bam = f"{RESULTS_DIR}/bowtie2_align/{{sample}}-{{genome}}.bam"
    output: 
        bam = f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam",
        flagstat_json = f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.flagstat.json"
    log: 
        f"{config['logdir']}/filter/{{sample}}-{{genome}}.log"
    params:
        exclude_flag = config['samtools_exclude_flag'],
        prefix = lambda wildcards: "-".join(wildcards)
    shell:
        """
        exec >> {log} 2>&1

        sambamba markdup -t {resources.threads} {input.bam} /dev/stdout | \
        samtools view -@ {resources.threads} -F 1804 -f 2 -q 30 -h -u /dev/stdin | \
        samtools sort -m 1536M -@ {resources.threads} -n /dev/stdin | \
        mbuffer -m 1G | \
        samtools fixmate -@ {resources.threads} -r /dev/stdin /dev/stdout | \
        samtools view -@ {resources.threads} -F 1804 -f 2 -h -u /dev/stdin | \
        samtools sort -m 1536M -@ {resources.threads} /dev/stdin | \
        mbuffer -m 1G | \
        samtools view -@ {resources.threads} -h | \
        grep -v "chrM" | \
        samtools view -@ {resources.threads} -b > {output.bam}

        samtools quickcheck -v {output.bam}
        samtools index -@ {resources.threads} {output.bam}
        samtools flagstat -@ {resources.threads} --output-fmt json {output.bam} > {output.flagstat_json}
        """

    # shell:
    #     """
    #     exec >> {log} 2>&1
    #     echo {params.prefix}
    #     sambamba markdup -t {resources.threads} {input.bam} /dev/stdout | \
    #     samtools view -@ {resources.threads} -F 1804 -f 2 -q 30 -h -u /dev/stdin | \
    #     samtools sort -@ {resources.threads} -n /dev/stdin | \
    #     samtools fixmate -@ {resources.threads} -r /dev/stdin /dev/stdout | \
    #     samtools view -@ {resources.threads} -F 1804 -f 2 -h -u /dev/stdin | \
    #     samtools sort -@ {resources.threads} /dev/stdin | \
    #     samtools view -@ {resources.threads} -h | grep -v "chrM" | samtools view -@ {resources.threads} -b > {output.bam}

    #     samtools quickcheck -v {output.bam}
    #     samtools index -@ {resources.threads} {output.bam}
    #     samtools flagstat -@ {resources.threads} --output-fmt json {output.bam} > {output.flagstat_json}
    #     """

