"""Z-score log10 bigWig signal per anchor, as the ENCODE4 cCRE pipeline does."""

import argparse
import sys

import numpy as np
import polars as pl

# bigtools bigwigaverageoverbed columns. mean0 counts uncovered bases as zero.
AVERAGES = ["name", "size", "bases", "sum", "mean0", "mean"]


def build_parser():
    """Build the zscore.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the zscore.py flags.
    """
    p = argparse.ArgumentParser(description="Z-score log10 bigWig signal per anchor.")
    p.add_argument("--output", required=True, help="Parquet to write.")
    return p


def zscores(signal):
    """Z-score log10 signal over the anchors that have signal.

    Parameters
    ----------
    signal : numpy.ndarray
        Mean signal per anchor, with uncovered bases counted as zero.

    Returns
    -------
    numpy.ndarray
        Z-score per anchor, from the mean and population standard deviation
        of log10 signal over anchors with signal. Anchors without signal
        get -10.
    """
    has_signal = signal > 0
    log_signal = np.log10(signal[has_signal])
    z = np.full(len(signal), -10.0)
    z[has_signal] = (log_signal - log_signal.mean()) / log_signal.std()
    return z


def main(argv=None):
    """Read bigwigaverageoverbed output from stdin and write z-scores by anchor.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    averages = pl.read_csv(
        sys.stdin.buffer,
        separator="\t",
        has_header=False,
        columns=[0, 4],
        new_columns=AVERAGES,
        schema_overrides={"name": pl.String, "mean0": pl.Float64},
    )
    signal = averages["mean0"].to_numpy()
    (
        pl.DataFrame(
            {"anchor": averages["name"], "zscore": zscores(signal), "signal": signal}
        )
        .with_columns(rank=pl.col("zscore").rank("min", descending=True))
        .sort("anchor")
        .write_parquet(args.output)
    )


if __name__ == "__main__":
    main()
