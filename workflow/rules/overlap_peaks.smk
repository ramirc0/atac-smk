rule overlap_peaks:
    input: 
        signal=get_macs3_pseudoreps_signal,
        pseudorep1="results/macs3_pseudoreps/{sample}_pseudorep1_peaks.narrowPeak", 
        pseudorep2="results/macs3_pseudoreps/{sample}_pseudorep2_peaks.narrowPeak",
    output: 
        overlap_peaks="results/overlap_peaks/{sample}_overlap_peaks.bed",
        signal=touch("results/overlap_peaks/{sample}.ok"),
    container: config['environments']['custom']
    log: "logs/overlap_peaks/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        bedtools intersect -a {input.pseudorep1} -b {input.pseudorep2} -u | sort -k1,1 -k2,2n > {output.overlap_peaks}
        """
