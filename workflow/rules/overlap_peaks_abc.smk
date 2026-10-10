rule overlap_peaks_abc: 
    input: 
        narrowpeak=rules.macs3_callpeak.output.narrowpeak,
        pseudorep1="results/macs3_pseudoreps/{sample}_pseudorep1_peaks.narrowPeak", 
        pseudorep2="results/macs3_pseudoreps/{sample}_pseudorep2_peaks.narrowPeak",
    output: 
        overlap_peaks="results/overlap_peaks_abc/{sample}.overlap_peaks.bed",
        signal=touch("results/overlap_peaks_abc/{sample}.ok"),
    software: SOFTWARE_ENV
    log: "logs/overlap_peaks_abc/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/overlap_peaks_abc.py {input.narrowpeak} {input.pseudorep1} {input.pseudorep2} {output.overlap_peaks}
        """
