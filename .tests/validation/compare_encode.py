"""Print the validation run's QC beside ENCODE's published values.

ENCODE ran older tools (MACS2 2.2.4, bowtie2 2.3.4.3, Picard) and stopped its
fragment length QC after 5M records, so values are close, not equal. The
integration tests compare step by step on identical inputs. Run from the
repository root after the validation run, which writes its QC table:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/validation/compare_encode.py
"""

import argparse
import subprocess
from pathlib import Path

import polars as pl
from zlib_ng import gzip_ng

ENCODE_DATA = Path("/data/projects/encode/data")
# The metrics of encode_qc.tsv, as qc_metrics.tsv names them after the rename.
COMPARED = [
    "nrf",
    "pbc1",
    "pbc2",
    "frac_mito",
    "frac_nfr",
    "nfr_mono_ratio",
    "frip_overlap",
    "frip_idr",
    "overlap_peaks",
    "idr_peaks",
]


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare QC with ENCODE's.")
    p.add_argument("--run", default="results/encode-validation", help="Run directory.")
    p.add_argument(
        "--encode", default=".tests/validation/encode_qc.tsv", help="ENCODE QC."
    )
    return p


def peaks_hit(a, b):
    """Return the fraction of peaks in ``a`` that overlap a peak in ``b``.

    Parameters
    ----------
    a, b : pathlib.Path
        Gzipped peak files.

    Returns
    -------
    float
        Overlapping peaks over all peaks of ``a``.
    """
    hits = subprocess.run(
        f"bedtools intersect -u -a {a} -b {b} | wc -l",
        shell=True,
        capture_output=True,
        text=True,
        check=True,
    )
    return int(hits.stdout) / peak_count(a)


def peak_count(path):
    """Count the peaks in a gzipped peak file.

    Parameters
    ----------
    path : pathlib.Path
        Gzipped peak file.

    Returns
    -------
    int
        Lines in the file.
    """
    with gzip_ng.open(path, "rt") as f:
        return sum(1 for _ in f)


def ours(run):
    """Read the run's QC table, keyed like ``encode_qc.tsv``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory. Its ``gather_qc`` rule wrote ``report/qc_metrics.tsv``.

    Returns
    -------
    dict of str to dict
        Per sample ID, the metrics ``encode_qc.tsv`` holds.
    """
    metrics = pl.read_csv(run / "report/qc_metrics.tsv", separator="\t").rename(
        {"frip": "frip_overlap"}
    )
    return {
        row["sample_id"]: {metric: row[metric] for metric in COMPARED}
        for row in metrics.iter_rows(named=True)
    }


def main(argv=None):
    """Print each metric per sample, ours beside ENCODE's, then peak agreement.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run)
    encode = pl.read_csv(args.encode, separator="\t", comment_prefix="#")
    measured = ours(run)
    for row in encode.iter_rows(named=True):
        sample = row["sample_id"]
        mine = measured[sample]
        print(f"\n{sample}\n{'metric':<16}{'ours':>12}{'ENCODE':>12}")
        for metric, value in mine.items():
            print(f"{metric:<16}{value:>12.4g}{row[metric]:>12.4g}")
        for peak_set in ("overlap", "idr"):
            path = run / f"{sample}/{peak_set}/{sample}.{peak_set}.bfilt.narrowPeak.gz"
            theirs = ENCODE_DATA / row[f"{peak_set}_file"]
            print(
                f"{peak_set} peaks: ours in ENCODE's {peaks_hit(path, theirs):.3f}, "
                f"ENCODE's in ours {peaks_hit(theirs, path):.3f}"
            )


if __name__ == "__main__":
    main()
