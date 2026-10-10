# Reads to a filtered BAM and its tagAlign.


rule fastp:
    input:
        r1=lambda w: R1_OF[w.sample],
        r2=lambda w: R2_OF[w.sample],
    output:
        r1=f"{OUTDIR}/{{sample}}/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/fastp/{{sample}}_R2.fastq.gz",
        html=f"{OUTDIR}/{{sample}}/fastp/{{sample}}.html",
        json=f"{OUTDIR}/{{sample}}/fastp/{{sample}}.json",
    log:
        f"{LOGDIR}/fastp/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fastp/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        fastp \
            -i {input.r1:q} \
            -I {input.r2:q} \
            -o {output.r1:q} \
            -O {output.r2:q} \
            -h {output.html:q} \
            -j {output.json:q} \
            -w {threads:q}
        """


# Run-independent: built once under results/refs/, reused across runs.
# conda's bowtie2-build lacks SAIS, so its size check alone picks a small .bt2
# index for GRCh38. --large-index keeps the .bt2l index the pipeline expects.
rule bowtie2_build:
    input:
        REFERENCE_FASTA,
    output:
        index=BOWTIE2_INDEX_FILES,
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
        chromsizes_bed=f"{REFDIR}/bowtie2_build/reference.chromsizes.bed",
    params:
        prefix=BOWTIE2_INDEX,
    log:
        f"{LOGDIR}/bowtie2_build.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_build.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools faidx -@ {threads:q} {input:q} -o - | cut -f 1,2 > {output.chromsizes:q}
        cat {output.chromsizes:q} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print $1,0,$2}}' \
            | sort -k1,1 -k2,2n \
            > {output.chromsizes_bed:q}
        bowtie2-build \
            --large-index \
            --threads {threads:q} \
            {input:q} \
            {params.prefix:q}
        """


rule bowtie2_align:
    input:
        index=BOWTIE2_INDEX_FILES,
        r1=f"{OUTDIR}/{{sample}}/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{sample}}/fastp/{{sample}}_R2.fastq.gz",
    output:
        bam=f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
        flagstat=f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.flagstat.json",
    params:
        k=lambda w: config["params"]["k"],
        index=BOWTIE2_INDEX,
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/bowtie2_align/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bowtie2_align/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bowtie2 \
            -X 2000 \
            --mm \
            -k {params.k:q} \
            --threads {threads:q} \
            -x {params.index:q} \
            -1 {input.r1:q} \
            -2 {input.r2:q} \
            | samtools view -@ {threads:q} -S -1 -b /dev/stdin \
            | samtools sort \
                -m 1536M \
                -@ {threads:q} \
                -T {params.tmp_dir:q}/{wildcards.sample:q}.$$ \
                -o {output.bam:q}
        samtools index -@ {threads:q} {output.bam:q}
        samtools flagstat -@ {threads:q} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


# -F excludes flag bits: 4 unmapped, 8 mate unmapped, 256 secondary, 512 failed QC,
# 1024 duplicate (1804 in total). -f 2 keeps proper pairs.
rule filter:
    input:
        f"{OUTDIR}/{{sample}}/bowtie2_align/{{sample}}.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
        flagstat=f"{OUTDIR}/{{sample}}/filter/{{sample}}.flagstat.json",
    params:
        exclude_flag=lambda w: config["params"]["samtools_exclude_flag"],
        include_flag=lambda w: config["params"]["samtools_include_flag"],
    log:
        f"{LOGDIR}/filter/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/filter/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        sambamba markdup -t {threads:q} {input:q} /dev/stdout \
            | samtools view -@ {threads:q} -F {params.exclude_flag:q} -f {params.include_flag:q} -q 30 -h -u /dev/stdin \
            | samtools sort -m 1536M -@ {threads:q} -n /dev/stdin \
            | mbuffer -m 1G \
            | samtools fixmate -@ {threads:q} -r /dev/stdin /dev/stdout \
            | samtools view -@ {threads:q} -F {params.exclude_flag:q} -f {params.include_flag:q} -h -u /dev/stdin \
            | samtools sort -m 1536M -@ {threads:q} /dev/stdin \
            | mbuffer -m 1G \
            | samtools view -@ {threads:q} -h \
            | grep -v "chrM" \
            | samtools view -@ {threads:q} -b \
            > {output.bam:q}

        samtools quickcheck -v {output.bam:q}
        samtools index -@ {threads:q} {output.bam:q}
        samtools flagstat -@ {threads:q} --output-fmt json {output.bam:q} > {output.flagstat:q}
        """


rule bam_to_tagalign:
    input:
        f"{OUTDIR}/{{sample}}/filter/{{sample}}.bam",
    output:
        bedpe=f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.bedpe",
        tagalign=f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign.gz",
    params:
        # The script gzips this file into the tagalign output.
        tagalign=f"{OUTDIR}/{{sample}}/bedpe_tagalign/{{sample}}.tagalign",
    log:
        f"{LOGDIR}/bam_to_tagalign/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/bam_to_tagalign/{{sample}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/bam_to_tagalign.py \
            --bam {input:q} \
            --bedpe {output.bedpe:q} \
            --tagalign {params.tagalign:q} \
            --threads {threads:q}
        """
