# MACS3 signal tracks as bigWig, as the ENCODE ATAC-seq pipeline makes them, and
# fold-change z-scores over the cCRE anchors.


# One job per track, so the three run in parallel. The bedGraphs are per million
# reads (--SPMR), so -S scales the p-value track back to read counts: tagAlign
# lines over a million, as ENCODE computes it. Fold change needs no scaling. The
# pileup track is not an ENCODE output.
rule signal:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        treat=f"{OUTDIR}/{{sample}}/peaks/{{sample}}_treat_pileup.bdg",
        control=f"{OUTDIR}/{{sample}}/peaks/{{sample}}_control_lambda.bdg",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        f"{OUTDIR}/{{sample}}/signal/{{sample}}.{{track}}.bigWig",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/signal/{{sample}}.{{track}}.txt",
    benchmark:
        f"{BENCHDIR}/signal/{{sample}}.{{track}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        case {wildcards.track:q} in
            fc)
                macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix track -m FE
                bedgraph="$tmp/track_FE.bdg"
                ;;
            pval)
                sval=$(python -m zlib_ng.gzip_ng -cd {input.tagalign:q} | wc -l | awk '{{ printf "%.6f", $1 / 1000000 }}')
                macs3 bdgcmp -t {input.treat:q} -c {input.control:q} --outdir "$tmp" --o-prefix track -m ppois -S "$sval"
                bedgraph="$tmp/track_ppois.bdg"
                ;;
            pileup)
                bedgraph={input.treat:q}
                ;;
        esac

        # ENCODE: clamp to the contigs, sort, and drop intervals that start inside
        # the previous one.
        bedtools slop -i "$bedgraph" -g {input.chromsizes:q} -b 0 \
            | bedClip stdin {input.chromsizes:q} /dev/stdout \
            | sort -k1,1 -k2,2n -T "$tmp" \
            | awk 'BEGIN {{ OFS = "\t" }} {{
                if (NR == 1 || prev_chr != $1 || prev_end <= $2) print
                prev_chr = $1; prev_end = $3
            }}' \
            > "$tmp/sorted.bdg"
        bedGraphToBigWig "$tmp/sorted.bdg" {input.chromsizes:q} {output:q}

        sync {output:q}
        """


# The ENCODE4 cCRE pipeline scores ATAC over unpadded anchors. SCREEN's ATAC
# z-scores come from ENCODE's fold-change track.
rule zscore:
    input:
        bigwig=f"{OUTDIR}/{{sample}}/signal/{{sample}}.fc.bigWig",
        anchors=ANCHORS,
    output:
        f"{OUTDIR}/{{sample}}/zscore/{{sample}}.parquet",
    log:
        f"{LOGDIR}/zscore/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/zscore/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bigwigaverageoverbed \
            --nthreads {threads:q} \
            {input.bigwig:q} \
            {input.anchors:q} \
            - \
            | python workflow/scripts/zscore.py --output {output:q}

        sync {output:q}
        """
