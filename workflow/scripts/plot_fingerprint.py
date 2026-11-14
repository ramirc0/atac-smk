"""Plot a fingerprint from deeptools plotFingerprint's raw bin counts.

Bins are sorted by read count. The curve is the fraction of reads in the
lowest-count fraction of bins. Uniform coverage follows the diagonal.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl

# The curve kept as data, at every 0.1% of bins. Smooth at any plot size.
CURVE_POINTS = 1001


def build_parser():
    """Build the plot_fingerprint.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_fingerprint.py flags.
    """
    p = argparse.ArgumentParser(description="Plot a fingerprint.")
    p.add_argument(
        "--counts", required=True, help="plotFingerprint --outRawCounts, one BAM."
    )
    p.add_argument("--title", default="", help="Figure title.")
    p.add_argument(
        "--curve",
        required=True,
        help="TSV to write: fraction_of_bins, fraction_of_reads.",
    )
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Draw the cumulative read fraction over the sorted bins, and keep it as data.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    raw = pl.read_csv(args.counts, separator="\t", comment_prefix="#", quote_char="'")
    counts = np.sort(raw[:, -1].to_numpy())
    bins = np.linspace(0, 1, counts.size)
    reads = np.cumsum(counts) / counts.sum()
    kept = np.linspace(0, 1, CURVE_POINTS)
    pl.DataFrame(
        {"fraction_of_bins": kept, "fraction_of_reads": np.interp(kept, bins, reads)}
    ).write_csv(args.curve, separator="\t")
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(bins, reads)
    ax.plot([0, 1], [0, 1], color=pk.PALETTES["plain"].foreground, linestyle=":")
    ax.set(xlabel="Fraction of bins", ylabel="Fraction of reads", title=args.title)
    pk.despine(ax)
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
