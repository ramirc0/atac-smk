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

rule chrombpnet_bias:
    input: 
        prep_nonpeaks_singal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        chromsizes = rules.bowtie2_build.output.chromsizes,
        peaks = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/peaks_no_blacklist.bed",
        negatives = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/negatives.bed",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    output:
        signal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_bias.continue",
        bias_model = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/models/bias.h5",
        bias_metrics = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/evaluation/bias_metrics.json"
    log: 
        f"{config['logdir']}/chrombpnet/bias/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json'],
        bias_threshold_factor = config['bias_threshold_factor']
    singularity: 'docker://clarity001/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1
        python -c 'import tensorflow as tf; print(f"Available compute devices: {{tf.config.list_physical_devices()}}")'
        OUTPUT_DIR=results/chrombpnet/{wildcards.sample}.{wildcards.genome}/bias_model/

        rm -rf $OUTPUT_DIR
        chrombpnet bias pipeline \
            -d ATAC \
            -itag {input.tagalign} \
            -g {params.fasta} \
            -c {input.chromsizes} \
            -p {input.peaks} \
            -n {input.negatives} \
            -fl {params.config_json} \
            -b {params.bias_threshold_factor} \
            -o $OUTPUT_DIR

        touch {output.signal}
        '''

rule chrombpnet_full:
    input: 
        prep_nonpeaks_singal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        chromsizes = rules.bowtie2_build.output.chromsizes,
        negatives = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/negatives.bed",
        peaks = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/peaks_no_blacklist.bed",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz",
        bias_model = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/models/bias.h5",
        bias_metrics = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/evaluation/bias_metrics.json"
    output:
        signal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_full.continue"
    log: 
        f"{config['logdir']}/chrombpnet/full/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json']
    singularity: 'docker://clarity001/chrombpnet:latest'
    shell:
        '''
        exec >> {log} 2>&1

        GT_THRESHOLD=$(python -c 'import json; bias_metrics = json.load(open("{input.bias_metrics}", "r")); gt_threshold = bias_metrics["counts_metrics"]["peaks"]["pearsonr"] > -0.5; print("Yes", end="") if gt_threshold else print("No", end="")')
        if [ "$GT_THRESHOLD" == "No" ]; then
            echo "The bias model has correlation below the threshold in peaks, skipping chrombpnet full."
            touch {output.signal}
            exit 0
        elif [ "$GT_THRESHOLD" == "Yes" ]; then
            echo "The bias model has correlation in peaks above the threshold, skipping chrombpnet full."
        fi

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
            -b {input.bias_model} \
            -o $OUTPUT_DIR

        touch {output.signal}
        '''
