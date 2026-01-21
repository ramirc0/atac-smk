# rule picard:
#     input: 
#         bam = f"{RESULTS_DIR}/filter/{{sample}}-{{genome}}.bam"
#     output: 
#         deduped_bam = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.bam", 
#         picard_report = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.txt",
#         flagstat_json = f"{RESULTS_DIR}/picard/{{sample}}-{{genome}}.flagstat.json"
#     params:
#         picard_tmpdir = config.get('picard_tmpdir', ".picard")
#     container: config["container"]
#     log: 
#         f"{config['logdir']}/picard/{{sample}}-{{genome}}.log"
#     shell:
#         """
#         exec >> {log} 2>&1
        
#         PICARD_TMP=$(mktemp -d -p {params.picard_tmpdir} picard_tmp.XXXXXX)
#         trap "rm -r $PICARD_TMP" EXIT

#         java -Dorg.xerial.snappy.tempdir=$PICARD_TMP -Djava.io.tmpdir=$PICARD_TMP -jar /usr/bin/picard.jar MarkDuplicates \
#             -INPUT {input.bam} \
#             -OUTPUT {output.deduped_bam} \
#             -METRICS_FILE {output.picard_report} \
#             -REMOVE_DUPLICATES TRUE \
#             -ASSUME_SORT_ORDER coordinate \
#             -USE_JDK_DEFLATER TRUE \
#             -USE_JDK_INFLATER TRUE \
#             -VALIDATION_STRINGENCY LENIENT
#         samtools quickcheck -v {output.deduped_bam}
#         samtools index -@ {resources.threads} {output.deduped_bam} {output.deduped_bam}.bai
#         samtools flagstat -@ --output-fmt json {resources.threads} {output.deduped_bam} > {output.flagstat_json}
#         """