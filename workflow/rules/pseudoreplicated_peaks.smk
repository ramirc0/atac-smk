# "Ideally we want Q = 0.05 for just the vanilla peak calls, and P = 0.01 for the pseudorep analysis"

rule pseudoreps:
    input: 
        bedpe = f"{RESULTS_DIR}/bedpe_tagalign/{{sample}}-{{genome}}.bedpe"
    output: 
        tagalign = [f"{RESULTS_DIR}/pseudoreps/{{sample}}-{{genome}}-pseudorep1.tagalign", f"{RESULTS_DIR}/pseudoreps/{{sample}}-{{genome}}-pseudorep2.tagalign"]
    container: config['container']
    log: f"{config['logdir']}/pseudoreps/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        ./workflow/rules/scripts/bedpe_to_pseudoreps.py {input.bedpe} {output.tagalign[0]} {output.tagalign[1]}
        """

rule macs3_pseudoreps:
    input:
        tagalign = f"{RESULTS_DIR}/pseudoreps/{{sample}}-{{genome}}-{{pseudorep}}.tagalign"
    output:
        narrowpeak = f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-{{pseudorep}}_peaks.narrowPeak"
    container: config['container']   
    log: f"{config['logdir']}/macs3_pseudoreps/{{sample}}-{{genome}}-{{pseudorep}}.log"
    shell:
        """
        exec >> {log} 2>&1
        macs3 callpeak -f BED -t {input} -n {wildcards.sample}-{wildcards.genome}-{wildcards.pseudorep} --outdir results/macs3_pseudoreps/ --nomodel --shift -75 --keep-dup all --extsize 150 -g hs --pvalue 0.01 -B
        """

rule idr:
    input: 
        narrowpeak = [f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep1_peaks.narrowPeak", f"results/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep2_peaks.narrowPeak"]
    output: 
        bed = f"{RESULTS_DIR}/idr/{{sample}}-{{genome}}.IDR.bed",
        png = f"{RESULTS_DIR}/idr/{{sample}}-{{genome}}.IDR.png"
    conda: config['conda_env']
    params:
        tmpdir = config['tmpdir']
    log: f"{config['logdir']}/idr/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        sort -k8,8nr {input.narrowpeak[0]} > {input.narrowpeak[0]}.sorted
        sort -k8,8nr {input.narrowpeak[1]} > {input.narrowpeak[1]}.sorted
        idr --samples {input.narrowpeak[0]}.sorted {input.narrowpeak[1]}.sorted --input-file-type narrowPeak --rank p.value --output-file /tmp/{wildcards.sample}-{wildcards.genome}.IDR.unsorted.bed --plot
        cat {params.tmpdir}/{wildcards.sample}-{wildcards.genome}.IDR.unsorted.bed | sort -k1,1 -k2,2n > {output.bed}
        mv {params.tmpdir}/{wildcards.sample}-{wildcards.genome}.IDR.unsorted.bed.png {output.png}
        """

rule overlap_peaks:
    input: 
        narrowpeak = [
            f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep1_peaks.narrowPeak", 
            f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep2_peaks.narrowPeak"
        ]
    output: 
        overlap_peaks = f"{RESULTS_DIR}/overlap_peaks/{{sample}}-{{genome}}.overlap_peaks.bed"
    container: config['container']
    log: f"{config['logdir']}/overlap_peaks/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        bedtools intersect -a {input.narrowpeak[0]} -b {input.narrowpeak[1]} -u | sort -k1,1 -k2,2n > {output.overlap_peaks}
        """

rule overlap_peaks_jill: 
    input: 
        narrowpeak = f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}.{{genome}}_peaks.narrowPeak",
        pseudorep1 = f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep1_peaks.narrowPeak", 
        pseudorep2 = f"{RESULTS_DIR}/macs3_pseudoreps/{{sample}}-{{genome}}-pseudorep2_peaks.narrowPeak",
    output: 
        overlap_peaks = f"{RESULTS_DIR}/overlap_peaks_jill/{{q}}/{{sample}}-{{genome}}.overlap_peaks.bed"
    container: config['container']
    log: f"{config['logdir']}/overlap_peaks_jill/{{q}}/{{sample}}-{{genome}}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/overlap_peaks_jill.py {input.narrowpeak} {input.pseudorep1} {input.pseudorep2} {output.overlap_peaks}
        """


