"""Compare each workflow step on S1 with ENCODE's own task script on the same inputs."""

import json
import subprocess

import polars as pl
import pytest
from conftest import (
    ANCHORS,
    BLACKLIST,
    CCRE_TOOLKIT,
    CHRSZ,
    FIXTURE,
    REFS,
    TSS,
    bigwig_stats,
    encode,
    first_difference,
    lines,
    run_rule,
)

PREFIX = "results-fixture/fixture/S1"


def encode_filter(run):
    """Filter S1's raw alignments with ENCODE's sambamba path.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.

    Returns
    -------
    pathlib.Path
        ENCODE output directory.
    """
    bam = run / "S1/bowtie2_align/S1.bam"
    # Its last step copies the input BAM into --out-dir, which must not hold the link.
    out = encode(
        "filter/S1",
        [bam, CHRSZ],
        f"python3 $(which encode_task_filter.py) {bam.name} --paired-end "
        "--multimapping 4 --dup-marker sambamba --mapq-thresh 30 "
        f"--filter-chrs chrM MT --chrsz {CHRSZ.name} --mito-chr-name chrM --nth 8 --mem-gb 8 --out-dir out",
    )
    return out / "out"


def encode_tagalign(run):
    """Convert S1's filtered BAM to a tagAlign with ENCODE's bam2ta.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.

    Returns
    -------
    pathlib.Path
        ENCODE's Tn5-shifted tagAlign.
    """
    bam = run / "S1/filter/S1.bam"
    out = encode(
        "bam2ta/S1",
        [bam],
        f"python3 $(which encode_task_bam2ta.py) {bam.name} --paired-end "
        "--mito-chr-name chrM --subsample 0 --nth 4",
    )
    return out / "S1.tn5.tagAlign.gz"


def alignments(bam):
    """List a BAM's alignments by name, flag, position, MAPQ and CIGAR.

    Parameters
    ----------
    bam : pathlib.Path
        BAM file.

    Returns
    -------
    list of str
        Sorted first six SAM fields of every record.
    """
    sam = subprocess.run(
        ["samtools", "view", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    return sorted("\t".join(r.split("\t")[:6]) for r in sam.splitlines())


def read_length(fastq):
    """Return ENCODE's read length: the longest of the first million reads.

    Parameters
    ----------
    fastq : pathlib.Path
        Gzipped FASTQ.

    Returns
    -------
    int
        Longest sequence among the first million reads.
    """
    return max(len(s) for s in lines(fastq)[1:4_000_000:4])


def nucleosomal_qc(path):
    """Parse a ``nucleosomal.qc`` file into its pass flags and values.

    Parameters
    ----------
    path : pathlib.Path
        File written in ENCODE's ``nucleosomal.qc`` format.

    Returns
    -------
    dict of str to tuple
        Per metric, the pass flag and, for ratio metrics, the value.
    """
    rows = [line.split("\t") for line in lines(path)]
    return {r[0]: (r[1], float(r[2])) if len(r) == 4 else (r[1],) for r in rows}


def test_filtered_bam_matches_encode(run):
    """Check that filtering keeps ENCODE's alignments.

    ENCODE's ``assign_multimappers.py`` never writes its last read-name
    group, so ours may hold one more read name than theirs.
    """
    ours = set(alignments(run / "S1/filter/S1.bam"))
    theirs = set(alignments(encode_filter(run) / "S1.nodup.no_chrM_MT.bam"))
    assert not theirs - ours, f"{len(theirs - ours)} of ENCODE's alignments missing"
    extra = {r.split("\t")[0] for r in ours - theirs}
    assert len(extra) <= 1, f"{len(extra)} extra read names"


def test_library_complexity_matches_encode(run):
    """Check NRF, PBC1 and PBC2 against ENCODE's on the same raw alignments."""
    qc = encode_filter(run) / "S1.dupmark.lib_complexity.qc"
    nrf, pbc1, pbc2 = map(float, lines(qc)[-1].split("\t")[4:])
    ours = json.loads((run / "S1/pbc_qc/S1.pbc.json").read_text())
    assert ours["NRF"] == pytest.approx(nrf, abs=1e-5)
    assert ours["PBC1"] == pytest.approx(pbc1, abs=1e-5)
    assert ours["PBC2"] == pytest.approx(pbc2, abs=1e-5)


def test_tagalign_matches_encode(run):
    """Check that the tagAlign has ENCODE's lines in ENCODE's order."""
    ours = lines(run / "S1/tagalign/S1.tagAlign.gz")
    assert first_difference(ours, lines(encode_tagalign(run))) is None


def test_pseudoreplicates_match_encode(run):
    """Check that pseudoreplicates keep mates together, shifted, split as ENCODE's."""
    ta = encode_tagalign(run)
    out = encode(
        "spr/S1",
        [ta],
        f"python3 $(which encode_task_spr.py) {ta.name} --paired-end "
        "--pseudoreplication-random-seed 0",
    )
    for pr in (1, 2):
        ours = lines(run / f"S1/pseudoreps/S1.pr{pr}.tagAlign.gz")
        theirs = lines(out / f"S1.tn5.pr{pr}.tagAlign.gz")
        assert first_difference(ours, theirs) is None


def test_pvalue_track_matches_encode(run):
    """Check the p-value bigWig against ENCODE's MACS2 track, within 2%."""
    ta = encode_tagalign(run)
    out = encode(
        "signal/S1",
        [ta, CHRSZ],
        f"python3 $(which encode_task_macs2_signal_track_atac.py) {ta.name} "
        f"--gensz hs --chrsz {CHRSZ.name} --pval-thresh 0.01 --smooth-win 150",
    )
    theirs = bigwig_stats(out / "S1.tn5.pval.signal.bigwig")
    ours = bigwig_stats(run / "S1/signal/S1.pval.bigWig")
    assert ours["mean"] == pytest.approx(theirs["mean"], rel=0.02)
    assert ours["max"] == pytest.approx(theirs["max"], rel=0.02)


def test_zscores_match_encode_ccre_pipeline(run):
    """Check fold-change z-scores against the ENCODE4 cCRE pipeline's script.

    UCSC's bigWigAverageOverBed rounds each anchor's mean to six significant
    digits before ENCODE z-scores it. bigtools does not, so values agree to
    1e-5, and anchors tied in ENCODE's ranking may split in ours.
    """
    bw = run / "S1/signal/S1.fc.bigWig"
    script = CCRE_TOOLKIT / "log-zscore-normalization.py"
    out = encode(
        "zscore/S1",
        [bw, ANCHORS, script],
        f"bigWigAverageOverBed {bw.name} {ANCHORS.name} averages.tab && "
        f"python3 {script.name} averages.tab > S1.zscore.txt",
    )
    # Each row holds an anchor, its z-score, signal and log10 signal.
    rows = sorted(line.split() for line in lines(out / "S1.zscore.txt"))
    z = pl.Series([float(r[1]) for r in rows])
    ours = pl.read_parquet(run / "S1/zscore/S1.parquet")
    assert first_difference(ours["anchor"].to_list(), [r[0] for r in rows]) is None
    assert ours["zscore"].to_list() == pytest.approx(z.to_list(), abs=1e-5)
    signal = [float(r[2]) for r in rows]
    assert ours["signal"].to_list() == pytest.approx(signal, rel=1e-5)
    # Read in our rank order, ENCODE's z-scores never rise.
    assert z.gather(ours["rank"].arg_sort()).diff().drop_nulls().le(0).all()


def peak_files(run):
    """Return S1's pseudoreplicate peaks, sample peaks and tagAlign.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.

    Returns
    -------
    list of pathlib.Path
        Pseudoreplicate 1 and 2 peaks, sample peaks, tagAlign.
    """
    return [
        run / "S1/peaks/S1.pr1.narrowPeak.gz",
        run / "S1/peaks/S1.pr2.narrowPeak.gz",
        run / "S1/peaks/S1.narrowPeak.gz",
        run / "S1/tagalign/S1.tagAlign.gz",
    ]


def encode_reproducibility(run, task):
    """Run ENCODE's IDR or naive overlap on S1's peaks.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    task : {"idr", "overlap"}
        ENCODE task.

    Returns
    -------
    pathlib.Path
        ENCODE output directory.
    """
    pr1, pr2, peaks, ta = peak_files(run)
    flags = "--idr-thresh 0.05 --idr-rank p.value" if task == "idr" else "--nonamecheck"
    return encode(
        f"{task}/S1",
        [pr1, pr2, peaks, ta, CHRSZ, BLACKLIST],
        f"python3 $(which encode_task_{task}.py) {pr1.name} {pr2.name} {peaks.name} "
        f"--prefix S1 --peak-type narrowPeak {flags} "
        f"--blacklist {BLACKLIST.name} --regex-bfilt-peak-chr-name 'chr[\\dXY]+' "
        f"--ta {ta.name} --chrsz {CHRSZ.name}",
    )


def frip(run, peak_set):
    """Return one FRiP row of S1's ``frip.tsv`` as reads in peaks over reads.

    Parameters
    ----------
    run : pathlib.Path
        Workflow run directory.
    peak_set : str
        Row name: ``peaks``, ``idr`` or ``overlap``.

    Returns
    -------
    float
        Fraction of tagAlign lines in the peak set.
    """
    rows = [line.split("\t") for line in lines(run / "S1/frip/S1.tsv")]
    row = next(r for r in rows if r[0] == peak_set)
    return int(row[1]) / int(row[2])


def test_peaks_match_encode(run):
    """Check MACS3 peak regions against ENCODE's MACS2 call on the same tagAlign.

    MACS3 3.0.5 fixed summit selection bugs (#741, #747, #748) that MACS2
    2.2.4 still has, so rows may differ in summit count and position.
    Regions may not.
    """
    ta = encode_tagalign(run)
    out = encode(
        "macs2/S1",
        [ta, CHRSZ],
        f"python3 $(which encode_task_macs2_atac.py) {ta.name} --gensz hs "
        f"--chrsz {CHRSZ.name} --pval-thresh 0.01 --smooth-win 150 "
        "--cap-num-peak 300000 --mem-gb 4",
    )
    ours = lines(run / "S1/peaks/S1.narrowPeak.gz")
    theirs = lines(out / "S1.tn5.pval0.01.300K.narrowPeak.gz")
    region = lambda line: tuple(line.split("\t")[:3])
    assert {region(x) for x in ours} == {region(x) for x in theirs}


def test_blacklist_filter_and_frip_match_encode(run):
    """Check the filtered sample peaks and their FRiP against ENCODE's."""
    _, _, peaks, ta = peak_files(run)
    out = encode(
        "post_call_peak/S1",
        [peaks, ta, CHRSZ, BLACKLIST],
        f"python3 $(which encode_task_post_call_peak_atac.py) {peaks.name} "
        f"--ta {ta.name} --peak-type narrowPeak --chrsz {CHRSZ.name} "
        f"--blacklist {BLACKLIST.name} --regex-bfilt-peak-chr-name 'chr[\\dXY]+'",
    )
    ours = lines(run / "S1/peaks/S1.bfilt.narrowPeak.gz")
    assert first_difference(ours, lines(out / "S1.bfilt.narrowPeak.gz")) is None
    theirs = float(lines(out / "S1.bfilt.frip.qc")[0])
    assert frip(run, "peaks") == pytest.approx(theirs, rel=1e-12)


def test_idr_matches_encode(run):
    """Check IDR peaks and their FRiP against ENCODE's on the same input peaks."""
    out = encode_reproducibility(run, "idr")
    ours = lines(run / "S1/idr/S1.idr.bfilt.narrowPeak.gz")
    theirs = lines(out / "S1.idr0.05.bfilt.narrowPeak.gz")
    assert first_difference(ours, theirs) is None
    theirs = float(lines(out / "S1.idr0.05.bfilt.frip.qc")[0])
    assert frip(run, "idr") == pytest.approx(theirs, rel=1e-12)


def test_overlap_matches_encode(run):
    """Check naive overlap peaks and their FRiP against ENCODE's."""
    out = encode_reproducibility(run, "overlap")
    ours = lines(run / "S1/overlap/S1.overlap.bfilt.narrowPeak.gz")
    theirs = lines(out / "S1.overlap.bfilt.narrowPeak.gz")
    assert first_difference(ours, theirs) is None
    theirs = float(lines(out / "S1.overlap.bfilt.frip.qc")[0])
    assert frip(run, "overlap") == pytest.approx(theirs, rel=1e-12)


def test_tss_enrichment_matches_encode(run):
    """Check the TSS enrichment score against ENCODE's metaseq code on the same BAM."""
    bam = run / "S1/filter/S1.bam"
    n = read_length(run / "S1/trim/S1_R1.fastq.gz")
    out = encode(
        "tss/S1",
        [bam, CHRSZ, TSS],
        f"python2 $(which encode_task_tss_enrich.py) --nodup-bam {bam.name} "
        f"--chrsz {CHRSZ.name} --tss {TSS.name} --read-len {n}",
    )
    theirs = float(lines(out / "S1.tss_enrich.qc")[0])
    ours = float(lines(run / "S1/tss_enrichment/S1.tss_enrich.qc")[0])
    assert ours == pytest.approx(theirs, rel=1e-9)


def test_fragment_length_qc_matches_encode(run):
    """Check the NFR fraction, NFR to mono-nucleosome ratio and peak checks."""
    bam = run / "S1/filter/S1.bam"
    out = encode(
        "fraglen/S1",
        [bam],
        f"python3 $(which encode_task_fraglen_stat_pe.py) --nodup-bam {bam.name} "
        "--picard-java-heap 4G",
    )
    ours = nucleosomal_qc(run / "S1/frag_len/S1.nucleosomal.qc")
    theirs = nucleosomal_qc(out / "S1.nucleosomal.qc")
    assert ours.keys() == theirs.keys()
    for metric, value in theirs.items():
        assert ours[metric] == pytest.approx(value, rel=1e-9)


def test_tagalign_job_fails_on_a_corrupt_bam(tmp_path):
    """Check that a corrupt BAM fails the tagAlign job instead of leaving empty files."""
    done = run_rule(
        tmp_path,
        "tagalign",
        f"{PREFIX}/tagalign/S1.tagAlign.gz",
        {f"{PREFIX}/filter/S1.bam": b"not a BAM"},
    )
    assert done.returncode != 0


def test_frip_job_fails_on_unreadable_peaks(tmp_path, run):
    """Check that a bedtools error fails the FRiP job instead of reporting 0."""
    files = {
        f"{PREFIX}/tagalign/S1.tagAlign.gz": run / "S1/tagalign/S1.tagAlign.gz",
        f"{PREFIX}/peaks/S1.bfilt.narrowPeak.gz": b"chr21\tnot\ta peak\n",
    }
    for path in ("idr/S1.idr", "overlap/S1.overlap"):
        files[f"{PREFIX}/{path}.bfilt.narrowPeak.gz"] = (
            run / f"S1/{path}.bfilt.narrowPeak.gz"
        )
    done = run_rule(tmp_path, "frip", f"{PREFIX}/frip/S1.tsv", files)
    assert done.returncode != 0


def test_idr_leaves_other_samples_scratch(tmp_path, run):
    """Check that S1's IDR job keeps S10's scratch files, whose names it prefixes."""
    scratch = tmp_path / ".tmp/S10_pseudorep1.sorted"
    files = {f"{PREFIX}/peaks/{p.name}": p for p in peak_files(run)[:3]}
    files["results-fixture/refs/bowtie2_build/reference.chromsizes"] = CHRSZ
    files[".tmp/S10_pseudorep1.sorted"] = b"S10 scratch\n"
    done = run_rule(tmp_path, "idr", f"{PREFIX}/idr/S1.idr.txt.gz", files)
    assert done.returncode == 0, done.stderr[-3000:]
    assert scratch.exists()


def qc_table(path):
    """Read a tab-separated ``name<TAB>value...`` QC file into a dict.

    Parameters
    ----------
    path : pathlib.Path
        QC file.

    Returns
    -------
    dict of str to list of str
        Fields after the name, per line.
    """
    return {r[0]: r[1:] for r in (line.split("\t") for line in lines(path))}


def test_annotation_enrichment_matches_encode(run):
    """Check read fractions in DHS, blacklist, promoter and enhancer regions."""
    ta = run / "S1/tagalign/S1.tagAlign.gz"
    regions = [FIXTURE / f"{a}.bed.gz" for a in ("ENCFF304XEX", "ENCFF140XLU")]
    regions += [FIXTURE / "ENCFF212UAV.bed.gz", BLACKLIST]
    dhs, prom, enh, bl = (p.name for p in regions)
    out = encode(
        "annot_enrich/S1",
        [ta, *regions],
        f"python3 $(which encode_task_annot_enrich.py) --ta {ta.name} --dnase {dhs} "
        f"--blacklist {bl} --prom {prom} --enh {enh}",
    )
    ours = qc_table(run / "S1/annot_enrich/S1.annot_enrich.qc")
    theirs = qc_table(out / "S1.annot_enrich.qc")
    assert ours.keys() == theirs.keys()
    for name, (reads, fraction) in theirs.items():
        assert ours[name][0] == reads
        assert float(ours[name][1]) == pytest.approx(float(fraction), abs=1e-6)


def test_mito_fraction_matches_encode(run):
    """Check frac_mito against ENCODE's bowtie2 and SAMstats counts."""
    bam = run / "S1/bowtie2_align/S1.bam"
    r1, r2 = (run / f"S1/trim/S1_R{i}.fastq.gz" for i in (1, 2))
    index = sorted((REFS.parent / "bowtie2_build_mito").glob("mito.*.bt2"))
    out = encode(
        "frac_mito/S1",
        [bam, r1, r2, CHRSZ, *index],
        f"python3 $(which encode_task_post_align.py) {r1.name} {bam.name} "
        f"--chrsz {CHRSZ.name} --mito-chr-name chrM --nth 8 --mem-gb 8 && "
        f"bowtie2 -k 5 -X2000 --mm --threads 8 -x mito -1 {r1.name} -2 {r2.name} "
        "| samtools sort -n -O sam - | SAMstats --sorted_sam_file - --outf mito.qc && "
        "python3 $(which encode_task_frac_mito.py) non_mito/*.samstats.qc mito.qc",
    )
    ours = qc_table(run / "S1/frac_mito/S1.frac_mito.qc")
    theirs = qc_table(next(out.glob("*.frac_mito.qc")))
    assert ours["non_mito_reads"] == theirs["non_mito_reads"]
    assert ours["mito_reads"] == theirs["mito_reads"]
    ours, theirs = (float(t["frac_mito_reads"][0]) for t in (ours, theirs))
    assert ours == pytest.approx(theirs, abs=1e-6)


def test_gc_bias_matches_encode(run):
    """Check Picard's normalized coverage per GC bin against ENCODE's run."""
    bam = run / "S1/filter/S1.bam"
    fasta = FIXTURE / "genome.fa"
    out = encode(
        "gc_bias/S1",
        [bam, fasta],
        f"python3 $(which encode_task_gc_bias.py) --nodup-bam {bam.name} "
        f"--ref-fa {fasta.name} --picard-java-heap 4G",
    )
    column = "NORMALIZED_COVERAGE"

    def coverage(path):
        rows = [x.split("\t") for x in lines(path) if x and not x.startswith("#")]
        i = rows[0].index(column)
        return [float(r[i]) for r in rows[1:]]

    ours = coverage(run / "S1/gc_bias/S1.gc.txt")
    assert ours == pytest.approx(coverage(out / "S1.gc.txt"), rel=1e-6)


def test_jsd_matches_encode(run):
    """Check the deeptools fingerprint metrics against ENCODE's, within 1%."""
    bam = run / "S1/filter/S1.bam"
    out = encode(
        "jsd/S1",
        [bam, BLACKLIST],
        f"python3 $(which encode_task_jsd.py) {bam.name} --blacklist {BLACKLIST.name} "
        "--mapq-thresh 30 --nth 8",
    )
    ours = [float(x) for x in lines(run / "S1/jsd/S1.jsd.qc")[1].split("\t")[1:]]
    theirs = [float(x) for x in lines(next(out.glob("*.jsd.qc")))[0].split("\t")]
    assert ours == pytest.approx(theirs, rel=0.01)
