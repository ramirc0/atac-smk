"""Plot Picard's GC bias metrics: normalized coverage and base quality by GC."""

import argparse
import io

import matplotlib.pyplot as plt
import poikilos as pk
import polars as pl


def build_parser():
    """Build the plot_gc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_gc.py flags.
    """
    p = argparse.ArgumentParser(description="Plot GC bias.")
    p.add_argument(
        "--metrics", required=True, help="Picard CollectGcBiasMetrics output."
    )
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Draw coverage and base quality by GC, as ENCODE's plot_gc does.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    # Picard puts a blank line between its header comments and the table.
    with open(args.metrics) as f:
        table = "".join(line for line in f if line.strip() and line[0] != "#")
    gc = pl.read_csv(io.StringIO(table), separator="\t")
    colors = pk.PALETTES["plain"].cycle
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(gc["GC"], gc["NORMALIZED_COVERAGE"], color=colors[0])
    ax.set(xlabel="GC (%)", ylabel="Normalized coverage", xlim=(0, 100))
    quality = ax.twinx()
    quality.plot(gc["GC"], gc["MEAN_BASE_QUALITY"], color=colors[1])
    quality.set_ylabel("Mean base quality", color=colors[1])
    windows = ax.twinx()
    windows.fill_between(
        gc["GC"], gc["WINDOWS"] / gc["WINDOWS"].sum(), color=colors[2], alpha=0.2
    )
    windows.set_axis_off()
    ax.set_title("GC bias (shaded: genome windows by GC)")
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
