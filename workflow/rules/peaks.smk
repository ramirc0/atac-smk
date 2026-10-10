# Pseudoreplicates, MACS3 peaks, IDR and naive overlap, as the ENCODE ATAC-seq
# pipeline does for one replicate.

import math

# ENCODE's MACS2 call for ATAC-seq: every Tn5-shifted mate is a tag, smoothed over
# smooth_win bp centred on its 5' end.
MACS3_ARGS = [
    "-f", "BED",
    "-g", "hs",
    "-p", str(config["params"]["pval_thresh"]),
    "--shift", str(-round(config["params"]["smooth_win"] / 2)),
    "--extsize", str(config["params"]["smooth_win"]),
    "--nomodel",
    "--keep-dup", "all",
    "--call-summits",
]

# ENCODE's peak file: p-value order, Peak_N names, cap_num_peak rows, clipped to
# the contigs. The cap sits in awk because head would close the pipe early.
FORMAT_PEAKS = r"""
        LC_ALL=C sort -k8,8gr "$tmp/peaks_peaks.narrowPeak" \
            | awk -v cap={params.cap:q} 'BEGIN {{ OFS = "\t" }}
                NR <= cap {{
                    $4 = "Peak_" NR
                    if ($2 < 0) $2 = 0
                    if ($3 < 0) $3 = 0
                    if ($10 == -1) $10 = $2 + int(($3 - $2 + 1) / 2.0)
                    print
                }}' \
            > "$tmp/capped"
        bedClip -truncate "$tmp/capped" {input.chromsizes:q} "$tmp/clipped"
        python -m zlib_ng.gzip_ng -cn "$tmp/clipped" > {output.peaks:q}
"""


# ENCODE's spr: whole pairs, shuffled with a keystream seeded by the uncompressed
# tagAlign size, then halved. The ENCODE image's OpenSSL 1.0.2 derives the key
# with MD5; -md md5 keeps that on newer OpenSSL.
rule pseudoreps:
    input:
        f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
    output:
        pr1=f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}.pr1.tagAlign.gz",
        pr2=f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}.pr2.tagAlign.gz",
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/pseudoreps/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/pseudoreps/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        python -m zlib_ng.gzip_ng -cd {input:q} > "$tmp/ta"
        lines=$(wc -l < "$tmp/ta")
        seed=$(wc -c < "$tmp/ta")
        sed 'N;s/\n/\t/' "$tmp/ta" \
            | shuf --random-source=<(openssl enc -aes-256-ctr -md md5 -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null) \
            | split -d -l $(( (lines / 2 + 1) / 2 )) - "$tmp/pr"
        unpair() {{
            awk 'BEGIN {{ OFS = "\t" }} {{ print $1, $2, $3, $4, $5, $6; print $7, $8, $9, $10, $11, $12 }}' "$1" \
                | python -m zlib_ng.gzip_ng -cn
        }}
        unpair "$tmp/pr00" > {output.pr1:q}
        unpair "$tmp/pr01" > {output.pr2:q}
        """


# ENCODE runs this call twice, for peaks and for signal. One call gives both.
rule macs3:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/tagalign/{{sample}}.tagAlign.gz",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        peaks=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.narrowPeak.gz",
        treat=temp(f"{OUTDIR}/{{sample}}/peaks/{{sample}}_treat_pileup.bdg"),
        control=temp(f"{OUTDIR}/{{sample}}/peaks/{{sample}}_control_lambda.bdg"),
    params:
        args=MACS3_ARGS,
        cap=lambda w: config["params"]["cap_num_peak"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/macs3/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        macs3 callpeak -t {input.tagalign:q} -n peaks --outdir "$tmp" {params.args:q} -B --SPMR
        mv "$tmp/peaks_treat_pileup.bdg" {output.treat:q}
        mv "$tmp/peaks_control_lambda.bdg" {output.control:q}
        """ + FORMAT_PEAKS


rule macs3_pseudorep:
    input:
        tagalign=f"{OUTDIR}/{{sample}}/pseudoreps/{{sample}}.{{pr}}.tagAlign.gz",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        peaks=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.{{pr}}.narrowPeak.gz",
    params:
        args=MACS3_ARGS,
        cap=lambda w: config["params"]["cap_num_peak"],
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/macs3_pseudorep/{{sample}}.{{pr}}.txt",
    benchmark:
        f"{BENCHDIR}/macs3_pseudorep/{{sample}}.{{pr}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        macs3 callpeak -t {input.tagalign:q} -n peaks --outdir "$tmp" {params.args:q}
        """ + FORMAT_PEAKS


# ENCODE's blacklist filter for every peak set: no contact with an excluded region,
# score capped at 1000, chr1-22, X and Y only.
rule blacklist_filter:
    input:
        peaks=f"{OUTDIR}/{{sample}}/{{path}}.narrowPeak.gz",
        blacklist=BLACKLIST,
    output:
        f"{OUTDIR}/{{sample}}/{{path}}.bfilt.narrowPeak.gz",
    log:
        f"{LOGDIR}/blacklist_filter/{{sample}}/{{path}}.txt",
    benchmark:
        f"{BENCHDIR}/blacklist_filter/{{sample}}/{{path}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bedtools intersect -nonamecheck -v -a {input.peaks:q} -b {input.blacklist:q} \
            | awk 'BEGIN {{ OFS = "\t" }} {{ if ($5 > 1000) $5 = 1000; print }}' \
            | awk '/chr[0-9XY]+([^A-Za-z0-9_]|$)/' \
            | python -m zlib_ng.gzip_ng -cn \
            > {output:q}
        """


# The sample's own peaks are the oracle list. --soft-idr-threshold reports every
# merged peak; the cut on column 12, global -log10 IDR, follows.
rule idr:
    input:
        pr1=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.pr1.narrowPeak.gz",
        pr2=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.pr2.narrowPeak.gz",
        peaks=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.narrowPeak.gz",
        chromsizes=f"{REFDIR}/bowtie2_build/reference.chromsizes",
    output:
        raw=f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.txt.gz",
        peaks=temp(f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.narrowPeak.gz"),
    params:
        threshold=lambda w: config["params"]["idr_thresh"],
        min_log10=lambda w: -math.log10(config["params"]["idr_thresh"]),
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/idr/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/idr/{{sample}}.tsv"
    software:
        IDR_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        idr \
            --samples {input.pr1:q} {input.pr2:q} \
            --peak-list {input.peaks:q} \
            --input-file-type narrowPeak \
            --output-file "$tmp/idr" \
            --rank p.value \
            --soft-idr-threshold {params.threshold:q} \
            --use-best-multisummit-IDR
        bedClip -truncate "$tmp/idr" {input.chromsizes:q} "$tmp/clipped"
        python -m zlib_ng.gzip_ng -cn "$tmp/clipped" > {output.raw:q}
        awk 'BEGIN {{ OFS = "\t" }} $12 >= {params.min_log10} {{
                if ($2 < 0) $2 = 0
                print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
            }}' "$tmp/clipped" \
            | sort \
            | uniq \
            | sort -grk8,8 \
            | cut -f 1-10 \
            | python -m zlib_ng.gzip_ng -cn \
            > {output.peaks:q}
        """


rule plot_idr:
    input:
        f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.txt.gz",
    output:
        svg=f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.svg",
        png=f"{OUTDIR}/{{sample}}/idr/{{sample}}.idr.png",
    params:
        threshold=lambda w: config["params"]["idr_thresh"],
    log:
        f"{LOGDIR}/plot_idr/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_idr/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python workflow/scripts/plot_idr.py \
            --input {input:q} \
            --threshold {params.threshold:q} \
            --title {wildcards.sample:q} \
            --output {output.svg:q}
        """


# ENCODE's naive overlap: sample peaks that overlap a peak in each pseudoreplicate,
# by at least half of either peak.
rule overlap:
    input:
        pr1=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.pr1.narrowPeak.gz",
        pr2=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.pr2.narrowPeak.gz",
        peaks=f"{OUTDIR}/{{sample}}/peaks/{{sample}}.narrowPeak.gz",
    output:
        temp(f"{OUTDIR}/{{sample}}/overlap/{{sample}}.overlap.narrowPeak.gz"),
    params:
        tmp_dir=TMP_DIR,
    log:
        f"{LOGDIR}/overlap/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/overlap/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d -p {params.tmp_dir:q})
        trap 'rm -rf "$tmp"' EXIT

        # ENCODE decompresses first; bedtools mishandles two gzipped inputs.
        python -m zlib_ng.gzip_ng -cd {input.pr1:q} > "$tmp/pr1"
        python -m zlib_ng.gzip_ng -cd {input.pr2:q} > "$tmp/pr2"
        python -m zlib_ng.gzip_ng -cd {input.peaks:q} > "$tmp/peaks"
        half() {{
            awk 'BEGIN {{ FS = OFS = "\t" }} {{ s1 = $3 - $2; s2 = $13 - $12; if ($21 / s1 >= 0.5 || $21 / s2 >= 0.5) print }}' \
                | cut -f 1-10 \
                | sort \
                | uniq
        }}
        bedtools intersect -nonamecheck -wo -a "$tmp/peaks" -b "$tmp/pr1" \
            | half \
            | bedtools intersect -nonamecheck -wo -a stdin -b "$tmp/pr2" \
            | half \
            | python -m zlib_ng.gzip_ng -cn \
            > {output:q}
        """
