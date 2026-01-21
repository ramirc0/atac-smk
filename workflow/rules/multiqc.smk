# rule multiqc:
#     input:
#         # macs3_output = expand(
#         #     f"{RESULTS_DIR}/macs3_callpeak/{{q}}/{{sample}}.{{genome}}_peaks.narrowPeak",
#         #     sample=SAMPLES,
#         #     genome=GENOMES,
#         #     q=QVALS
#         # ),
#         chrombpnet_output = expand(
#             [
#                 f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_bias.continue",
#                 f"{RESULTS_DIR}/chrombpnet/{{sample}}.{{genome}}/chrombpnet_prep_nonpeaks.continue"
#             ],
#             sample=SAMPLES,
#             genome=GENOMES
#         )
#     output:
#         datadir = directory(f"{RESULTS_DIR}/multiqc/multiqc_data"),
#         html_report = f"{RESULTS_DIR}/multiqc/multiqc_report.html",
#     params: 
#         multiqc_conf = config["input_files"]['multiqc_config'],
#         macs3_qvalue = config['qvals'][0]
#     container: config['container']
#     log: f"{config['logdir']}/mutliqc.log"
#     shell:
#         """
#         exec >> {log} 2>&1
#         echo "$(date --iso=minutes): Started multiqc"

#         python workflow/rules/scripts/run_multiqc.py --multiqc_config {params.multiqc_conf} --outdir {RESULTS_DIR} --macs3_qvalue {params.macs3_qvalue}

#         echo "$(date --iso=minutes): Finished multiqc"
#         """
