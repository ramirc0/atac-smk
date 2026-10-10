"""Print the integration run's QC beside ENCODE's published values.

ENCODE ran older tools (MACS2 2.2.4, bowtie2 2.3.4.3, Picard) and stopped its
fragment length QC after 5M records, so values are close, not equal. The
fixture tests compare step by step on identical inputs. Run from the
repository root after the integration run:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/compare_encode.py
"""

import argparse
import json
import subprocess
from pathlib import Path

import polars as pl
from zlib_ng import gzip_ng

ENCODE_DATA = Path("/data/projects/encode/data")


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare QC with ENCODE's.")
    p.add_argument("--run", default="results/encode-integration", help="Run directory.")
    p.add_argument(
        "--encode", default=".tests/integration/encode_qc.tsv", help="ENCODE QC."
    )
    return p


def table(path):
    """Read a ``name<TAB>value...`` QC file into a dict of its first values.

    Parameters
    ----------
    path : pathlib.Path
        QC file.

    Returns
    -------
    dict of str to str
        First field after the name, per line.
    """
    return dict(line.split("\t")[:2] for line in path.read_text().splitlines())


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


def ours(run, sample):
    """Collect one sample's QC from the run, keyed like ``encode_qc.tsv``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.

    Returns
    -------
    dict of str to float
        QC values.
    """
    d = run / sample
    pbc = json.loads((d / f"pbc_qc/{sample}.pbc.json").read_text())
    nucleosomal = {
        k: v.split("\t")[1]
        for k, v in table(d / f"frag_len/{sample}.nucleosomal.qc").items()
    }
    frip = pl.read_csv(d / f"frip/{sample}.tsv", separator="\t")
    frip = dict(zip(frip["peak_set"], frip["reads_in_peaks"] / frip["reads"]))
    return {
        "nrf": pbc["NRF"],
        "pbc1": pbc["PBC1"],
        "pbc2": pbc["PBC2"],
        "frac_mito": float(
            table(d / f"frac_mito/{sample}.frac_mito.qc")["frac_mito_reads"]
        ),
        "frac_nfr": float(nucleosomal["Fraction of reads in NFR"]),
        "nfr_mono_ratio": float(nucleosomal["NFR / mono-nuc reads"]),
        "frip_overlap": frip["overlap"],
        "frip_idr": frip["idr"],
        "overlap_peaks": peak_count(
            d / f"overlap/{sample}.overlap.bfilt.narrowPeak.gz"
        ),
        "idr_peaks": peak_count(d / f"idr/{sample}.idr.bfilt.narrowPeak.gz"),
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
    for row in encode.iter_rows(named=True):
        sample = row["sample_id"]
        mine = ours(run, sample)
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
