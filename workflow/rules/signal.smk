# MACS3 signal tracks as bigWig, as the ENCODE ATAC-seq pipeline makes them.


# The bedGraphs are per million reads (--SPMR), so -S scales the p-value track back
# to read counts: tagAlign lines over a million, as ENCODE computes it. Fold change
# needs no scaling. The pileup bigWig is not an ENCODE output.
rule signal:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        treat=f"{OUTDIR}/{{sample}}/peaks/{{sample}}_treat_pileup.bdg",
        control=f"{OUTDIR}/{{sample}}/peaks/{{sample}}_control_lambda.bdg",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        fc=f"{OUTDIR}/{{sample}}/signal/{{sample}}.fc.bigWig",
        pval=f"{OUTDIR}/{{sample}}/signal/{{sample}}.pval.bigWig",
        pileup=f"{OUTDIR}/{{sample}}/signal/{{sample}}.pileup.bigWig",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/signal/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/signal/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        # ENCODE: clamp to the contigs, sort, and drop intervals that start inside
        # the previous one.
        bigwig() {{
            bedtools slop -i "$1" -g {input.chromsizes:q} -b 0 \
                | bedClip stdin {input.chromsizes:q} /dev/stdout \
                | sort -k1,1 -k2,2n -T "$tmp" \
                | awk 'BEGIN {{ OFS = "\t" }} {{
                    if (NR == 1 || prev_chr != $1 || prev_end <= $2) print
                    prev_chr = $1; prev_end = $3
                }}' \
                > "$tmp/sorted.bdg"
            bedGraphToBigWig "$tmp/sorted.bdg" {input.chromsizes:q} "$2"
        }}

        sval=$(python -m zlib_ng.gzip_ng -cd {input.tagalign:q} | wc -l | awk '{{ printf "%.6f", $1 / 1000000 }}')
        macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix sample -m FE
        macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix sample -m ppois -S "$sval"
        bigwig "$tmp/sample_FE.bdg" {output.fc:q}
        bigwig "$tmp/sample_ppois.bdg" {output.pval:q}
        bigwig {input.treat:q} {output.pileup:q}
        """
