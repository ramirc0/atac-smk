"""Write each sample's report figures, interactive where the data are kept."""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from qc_report import (
    gc_bias_figure,
    idr_signal_figure,
    read_picard_metrics,
    read_quality_figure,
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


def png_href(run, sample, step, report_dir):
    """Link a step's PNG for one sample from the report's directory.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.
    step : str
        Step directory, such as ``tss_enrichment``.
    report_dir : pathlib.Path
        Directory of the report HTML.

    Returns
    -------
    str
        Relative path with forward slashes.
    """
    png = run / sample / step / f"{sample}.png"
    return Path(os.path.relpath(png, report_dir)).as_posix()


def write_panels(sample, run, report_dir, idr_threshold):
    """Write one sample's ``samples/<sample>.js``.

    The TSS profile, fragment lengths and fingerprint are not kept as data,
    so their panels show the pipeline's PNG.

    Parameters
    ----------
    sample : str
        Sample ID.
    run : pathlib.Path
        Run directory.
    report_dir : pathlib.Path
        Directory of the report HTML. Image paths are relative to it.
    idr_threshold : float
        IDR threshold, for the IDR panel's legend.
    """
    d = run / sample
    write_sample_panels(
        report_dir / "samples" / f"{sample}.js",
        sample,
        [
            ("Read quality", read_quality_figure(d / f"fastp/{sample}.json")),
            ("TSS enrichment", png_href(run, sample, "tss_enrichment", report_dir)),
            ("Fragment lengths", png_href(run, sample, "frag_len", report_dir)),
            (
                "GC bias",
                gc_bias_figure(read_picard_metrics(d / f"gc_bias/{sample}.gc.txt")),
            ),
            ("Fingerprint", png_href(run, sample, "jsd", report_dir)),
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
