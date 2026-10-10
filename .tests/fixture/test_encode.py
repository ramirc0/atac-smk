"""Compare each workflow step on S1 with ENCODE's own task script on the same inputs."""

import json
import subprocess

import pytest
from conftest import (
    BLACKLIST,
    CHRSZ,
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
    ours = lines(run / "S1/bedpe_tagalign/S1.tagalign.gz")
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
        ours = lines(run / f"S1/pseudoreps/S1_pseudorep{pr}.tagalign")
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
    ours = bigwig_stats(run / "S1/macs3_pvalue_signal/S1_ppois.bw")
    assert ours["mean"] == pytest.approx(theirs["mean"], rel=0.02)
    assert ours["max"] == pytest.approx(theirs["max"], rel=0.02)


def test_idr_matches_encode(run):
    """Check IDR peaks against ENCODE's on the same pseudoreplicate and sample peaks."""
    pr1, pr2 = (
        run / f"S1/macs3_pseudoreps/S1_pseudorep{i}_peaks.narrowPeak" for i in (1, 2)
    )
    peaks = run / "S1/macs3_callpeak/S1_peaks.narrowPeak"
    ta = run / "S1/bedpe_tagalign/S1.tagalign.gz"
    out = encode(
        "idr/S1",
        [pr1, pr2, peaks, ta, CHRSZ, BLACKLIST],
        f"python3 $(which encode_task_idr.py) {pr1.name} {pr2.name} {peaks.name} "
        "--prefix S1 --peak-type narrowPeak --idr-thresh 0.05 --idr-rank p.value "
        f"--blacklist {BLACKLIST.name} --regex-bfilt-peak-chr-name 'chr[\\dXY]+' "
        f"--ta {ta.name} --chrsz {CHRSZ.name}",
    )
    ours = lines(run / "S1/idr/S1.IDR.bed")
    theirs = lines(out / "S1.idr0.05.bfilt.narrowPeak.gz")
    assert first_difference(ours, theirs) is None


def test_tss_enrichment_matches_encode(run):
    """Check the TSS enrichment score against ENCODE's metaseq code on the same BAM."""
    bam = run / "S1/filter/S1.bam"
    n = read_length(run / "S1/fastp/S1_R1.fastq.gz")
    out = encode(
        "tss/S1",
        [bam, CHRSZ, TSS],
        f"python2 $(which encode_task_tss_enrich.py) --nodup-bam {bam.name} "
        f"--chrsz {CHRSZ.name} --tss {TSS.name} --read-len {n}",
    )
    theirs = float(lines(out / "S1.tss_enrich.qc")[0])
    ours = json.loads((run / "S1/tss_enrichment/S1.txt").read_text())
    assert ours["tss_enrichment_score"] == pytest.approx(theirs, rel=1e-9)


def test_fragment_length_qc_matches_encode(run):
    """Check the NFR fraction, NFR to mono-nucleosome ratio and peak checks."""
    bam = run / "S1/filter/S1.bam"
    out = encode(
        "fraglen/S1",
        [bam],
        f"python3 $(which encode_task_fraglen_stat_pe.py) --nodup-bam {bam.name} "
        "--picard-java-heap 4G",
    )
    ours = nucleosomal_qc(run / "S1/frag_len/S1.txt")
    theirs = nucleosomal_qc(out / "S1.nucleosomal.qc")
    assert ours.keys() == theirs.keys()
    for metric, value in theirs.items():
        assert ours[metric] == pytest.approx(value, rel=1e-9)


def test_tagalign_job_fails_on_a_corrupt_bam(tmp_path):
    """Check that a corrupt BAM fails the tagAlign job instead of leaving empty files."""
    done = run_rule(
        tmp_path,
        "bam_to_tagalign",
        f"{PREFIX}/bedpe_tagalign/S1.tagalign.gz",
        {f"{PREFIX}/filter/S1.bam": b"not a BAM"},
    )
    assert done.returncode != 0


def test_frip_job_fails_on_unreadable_peaks(tmp_path, run):
    """Check that a bedtools error fails the FRiP job instead of reporting 0."""
    done = run_rule(
        tmp_path,
        "frip_all",
        f"{PREFIX}/frip_all/S1.txt",
        {
            f"{PREFIX}/bedpe_tagalign/S1.tagalign.gz": (
                run / "S1/bedpe_tagalign/S1.tagalign.gz"
            ),
            f"{PREFIX}/macs3_callpeak/S1_peaks.narrowPeak": b"chr21\tnot\ta peak\n",
        },
    )
    assert done.returncode != 0


def test_idr_leaves_other_samples_scratch(tmp_path, run):
    """Check that S1's IDR job keeps S10's scratch files, whose names it prefixes."""
    scratch = tmp_path / ".tmp/S10_pseudorep1.sorted"
    files = {
        f"{PREFIX}/macs3_pseudoreps/S1_pseudorep{i}_peaks.narrowPeak": (
            run / f"S1/macs3_pseudoreps/S1_pseudorep{i}_peaks.narrowPeak"
        )
        for i in (1, 2)
    }
    files[".tmp/S10_pseudorep1.sorted"] = b"S10 scratch\n"
    done = run_rule(tmp_path, "idr", f"{PREFIX}/idr/S1.IDR.bed", files)
    assert done.returncode == 0, done.stderr[-3000:]
    assert scratch.exists()
