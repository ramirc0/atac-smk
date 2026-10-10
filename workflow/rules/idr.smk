rule idr:
    input: 
        signal=get_macs3_pseudoreps_signal,
        pseudorep1="results/macs3_pseudoreps/{sample}_pseudorep1_peaks.narrowPeak", 
        pseudorep2="results/macs3_pseudoreps/{sample}_pseudorep2_peaks.narrowPeak",
    output: 
        bed="results/idr/{sample}.IDR.bed",
        png="results/idr/{sample}.IDR.png",
        signal=touch("results/idr/{sample}.ok"),
    software: IDR_ENV
    params:
        tmp_dir=config['input']['tmp_dir']
    log: "logs/idr/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1

        sort -k8,8nr {input.pseudorep1} > {params.tmp_dir}/{wildcards.sample}_pseudorep1.sorted
        sort -k8,8nr {input.pseudorep2} > {params.tmp_dir}/{wildcards.sample}_pseudorep2.sorted
        idr --samples {params.tmp_dir}/{wildcards.sample}_pseudorep1.sorted \
                      {params.tmp_dir}/{wildcards.sample}_pseudorep2.sorted \
            --input-file-type narrowPeak \
            --rank p.value \
            --output-file {params.tmp_dir}/{wildcards.sample}.IDR.unsorted.bed \
            --plot
        cat {params.tmp_dir}/{wildcards.sample}.IDR.unsorted.bed | sort -k1,1 -k2,2n > {output.bed}
        mv {params.tmp_dir}/{wildcards.sample}.IDR.unsorted.bed.png {output.png}
        rm {params.tmp_dir}/{wildcards.sample}*
        """
