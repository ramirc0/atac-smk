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
    singularity: config['chrombpnet_container']
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
    singularity: config['chrombpnet_container']
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

rule choose_best_bias_model:
    input:
        bias_models = expand(
            f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/models/bias.h5",
            sample=SAMPLES,
            genome=GENOMES
        ),
        bias_metrics = expand(
            f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/bias_model/evaluation/bias_metrics.json",
            sample=SAMPLES,
            genome=GENOMES
        )
    output:
        best_bias_model = f"{RESULTS_DIR}/chrombpnet/best_bias_model.h5",
        best_bias_metrics = f"{RESULTS_DIR}/chrombpnet/best_bias_metrics.json"
    log: f"{config['logdir']}/chrombpnet/best_bias_model.log"
    singularity: config['chrombpnet_container']
    run:
        import sys
        import json
        import shutil
        logfile = open(log[0], "w")
        sys.stdout = logfile
        sys.stdout = logfile

        best_pearsonr = None
        best_model_path = None
        best_metrics_path = None

        for bias_model, bias_metrics_path in zip(input.bias_models, input.bias_metrics):
            with open(bias_metrics_path, "r") as f:
                bias_metrics = json.load(f)
                peaks_pearsonr = bias_metrics["counts_metrics"]["peaks"]["pearsonr"]
                if best_pearsonr is None or peaks_pearsonr > best_pearsonr:
                    print("The new best peaks pearsonr is {}".format(peaks_pearsonr))
                    best_pearsonr = peaks_pearsonr
                    best_model_path = bias_model
                    best_metrics_path = bias_metrics_path

        print("The best bias model is {} with a peaks pearsonr of {}".format(best_model_path, best_pearsonr))
        shutil.copy(best_model_path, output.best_bias_model)
        shutil.copy(best_metrics_path, output.best_bias_metrics)
        print("Copied best bias model to {}".format(output.best_bias_model))
        print("Copied best bias metrics to {}".format(output.best_bias_metrics))


rule chrombpnet_full:
    input: 
        prep_nonpeaks_singal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue",
        chromsizes = rules.bowtie2_build.output.chromsizes,
        negatives = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/negatives.bed",
        peaks = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/peaks_no_blacklist.bed",
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz",
        bias_model = f"{RESULTS_DIR}/chrombpnet/best_bias_model.h5"
    output:
        signal = f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_full.continue"
    log: 
        f"{config['logdir']}/chrombpnet/full/{{sample}}.{{genome}}.log"
    params:
        fasta = lambda wildcards: config['genomes'].get(wildcards.genome),
        config_json = config["input_files"]['chrombpnet_config_json']
    singularity: config['chrombpnet_container']
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
            -b {input.bias_model} \
            -o $OUTPUT_DIR

        touch {output.signal}
        '''