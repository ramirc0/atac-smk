"""Write each sample's report figures, from the data each QC step keeps."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import polars as pl
from fragment_ranges import MONO_LOWER, MONO_UPPER, NFR_UPPER
from qc_report import (
    fingerprint_figure,
    gc_bias_figure,
    idr_signal_figure,
    insert_size_figure,
    read_picard_metrics,
    read_quality_figure,
    tss_profile_figure,
    write_sample_panels,
)


def build_parser():
    """Build the qc_figures.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_figures.py flags.
    """
    p = argparse.ArgumentParser(description="Write each sample's report figures.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--report-dir", required=True, help="Directory of the report HTML.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--idr-threshold", type=float, default=0.05, help="IDR threshold.")
    p.add_argument("--threads", type=int, default=1, help="Samples drawn at once.")
    return p


def write_panels(sample, run, report_dir, idr_threshold):
    """Write one sample's ``samples/<sample>.js``.

    Parameters
    ----------
    sample : str
        Sample ID.
    run : pathlib.Path
        Run directory.
    report_dir : pathlib.Path
        Directory of the report HTML.
    idr_threshold : float
        IDR threshold, for the IDR panel's legend.
    """
    d = run / sample
    histogram = read_picard_metrics(
        d / f"frag_len/{sample}.insert_size.txt", histogram=True
    )
    write_sample_panels(
        report_dir / "samples" / f"{sample}.js",
        sample,
        [
            ("Read quality", read_quality_figure(d / f"fastp/{sample}.json")),
            (
                "TSS enrichment",
                tss_profile_figure(
                    pl.read_csv(
                        d / f"tss_enrichment/{sample}.tss_profile.tsv", separator="\t"
                    )
                ),
            ),
            (
                "Fragment lengths",
                insert_size_figure(
                    histogram,
                    bands=[("NFR", 0, NFR_UPPER), ("Mono-nuc", MONO_LOWER, MONO_UPPER)],
                ),
            ),
            (
                "GC bias",
                gc_bias_figure(read_picard_metrics(d / f"gc_bias/{sample}.gc.txt")),
            ),
            (
                "Fingerprint",
                fingerprint_figure(
                    pl.read_csv(d / f"jsd/{sample}.fingerprint.tsv", separator="\t")
                ),
            ),
            (
                "IDR",
                idr_signal_figure(d / f"idr/{sample}.idr.txt.gz", idr_threshold),
            ),
        ],
    )


def main(argv=None):
    """Write every sample's figures.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    draw = partial(
        write_panels,
        run=Path(args.run_dir),
        report_dir=Path(args.report_dir),
        idr_threshold=args.idr_threshold,
    )
    with ProcessPoolExecutor(args.threads) as pool:
        list(pool.map(draw, args.samples))


if __name__ == "__main__":
    main()
