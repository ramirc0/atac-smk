rule macs3_callpeak:
    input:
        tagalign = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.tagalign.gz"
    output:
        narrowpeak = f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}.{{genome}}_peaks.narrowPeak"
    container: config["container"]
    log: 
        f"{config['logdir']}/macs3_callpeak/{{q}}/{{sample}}.{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input.tagalign} -n {wildcards.sample}.{wildcards.genome} --outdir {RESULTS_DIR}/macs3_callpeak/{wildcards.q}/ --nomodel --shift -75 --keep-dup all --extsize 150 -g hs -q {wildcards.q}
        """

