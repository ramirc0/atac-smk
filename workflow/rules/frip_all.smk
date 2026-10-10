rule frip_all:
    input: 
        peaks_signal=rules.macs3_callpeak.output.signal,
        tagalign_signal=rules.bam_to_tagalign.output.signal,
        peaks=rules.macs3_callpeak.output.narrowpeak,
        tagalign=rules.bam_to_tagalign.output.tagalign,
    output:
        frip_report = "results/frip_all/{sample}.txt",
        signal=touch("results/frip_all/{sample}.ok"),
    software: SOFTWARE_ENV
    log: 
        "logs/frip_all/{sample}.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/scripts/calculate_frip.py --tagalign {input.tagalign} --peaks {input.peaks} --output {output.frip_report}
        """ 
