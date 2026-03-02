rule chrombpnet_prep_nonpeaks:
    input: 
        peaks = f"{RESULTS_DIR}/macs3_callpeak/{MIN_QVAL}/{{sample}}.{{genome}}_peaks.narrowPeak",
        chromsizes = rules.bowtie2_build.output.chromsizes,
    output:
        signal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        negatives = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/negatives.bed",
        peaks = temp(f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/peaks_no_blacklist.bed")
    log: 
        f"{config['logdir']}/chrombpnet/prep_nonpeaks/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json'],
        blacklist = config["input_files"]['blacklist']
    singularity: 'docker://clarity001/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1
        python -c 'import tensorflow as tf; print(f"Available compute devices: {{tf.config.list_physical_devices()}}")'

        bedtools slop -i {params.blacklist} -g {input.chromsizes} -b 1057 | \
            bedtools intersect -v -a {input.peaks} -b - > {RESULTS_DIR}/chrombpnet/{wildcards.sample}.{wildcards.genome}/peaks_no_blacklist.bed

        chrombpnet prep nonpeaks \
            -g {params.fasta} \
            -p {RESULTS_DIR}/chrombpnet/{wildcards.sample}.{wildcards.genome}/peaks_no_blacklist.bed \
            -c {input.chromsizes} \
            -fl {params.config_json} \
            -br {params.blacklist} \
            -o results/chrombpnet/{wildcards.sample}.{wildcards.genome}
        
        mv results/chrombpnet/{wildcards.sample}.{wildcards.genome}_negatives.bed results/chrombpnet/{wildcards.sample}.{wildcards.genome}/negatives.bed
        mv results/chrombpnet/{wildcards.sample}.{wildcards.genome}_auxiliary/ results/chrombpnet/{wildcards.sample}.{wildcards.genome}/prep_nonpeaks_auxiliary/
        
        touch {output.signal}
        '''

rule chrombpnet_full:
    input:
        prep_nonpeaks_singal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        chromsizes = rules.bowtie2_build.output.chromsizes,
        negatives = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/negatives.bed",
        peaks = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/peaks_no_blacklist.bed",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz",
    output:
        signal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_full.continue"
    log: 
        f"{config['logdir']}/chrombpnet/full/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        chrombpnet_bias_model = config['input_files']['chrombpnet_bias_model'],
        config_json = config["input_files"]['chrombpnet_config_json']
    singularity: 'docker://clarity001/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1

        python -c 'import tensorflow as tf; print(f"Available compute devices: {{tf.config.list_physical_devices()}}")'
        OUTPUT_DIR=results/chrombpnet/{wildcards.sample}.{wildcards.genome}/full/

        rm -rf $OUTPUT_DIR
        chrombpnet pipeline \
            -d ATAC \
            -itag {input.tagalign} \
            -g {params.fasta} \
            -c {input.chromsizes} \
            -p {input.peaks} \
            -n {input.negatives} \
            -fl {params.config_json} \
            -b {params.chrombpnet_bias_model} \
            -o $OUTPUT_DIR

        touch {output.signal}
        '''
