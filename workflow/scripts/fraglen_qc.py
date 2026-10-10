"""Score and plot fragment lengths from Picard's insert size histogram, as ENCODE does."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
from scipy.signal import find_peaks_cwt

# ENCODE's ranges in bp: NFR below 150, mono-nucleosome above 150 up to 300.
NFR_UPPER = 150
MONO_LOWER = 150
MONO_UPPER = 300
PEAK_RANGES = [
    ("Presence of NFR peak", 20, 90),
    ("Presence of Mono-Nuc peak", 120, 250),
    ("Presence of Di-Nuc peak", 300, 500),
]


def build_parser():
    """Build the fraglen_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the fraglen_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Score and plot fragment lengths.")
    p.add_argument(
        "--histogram", required=True, help="Picard CollectInsertSizeMetrics output."
    )
    p.add_argument("--qc", required=True, help="nucleosomal.qc to write.")
    p.add_argument(
        "--plot", required=True, help="Distribution SVG to write (PNG beside it)."
    )
    return p


def read_histogram(path):
    """Read the histogram block of a Picard insert size metrics file.

    Parameters
    ----------
    path : str
        Picard CollectInsertSizeMetrics output.

    Returns
    -------
    numpy.ndarray
        One row per insert size: the size, then the count per pair
        orientation. The first count column is ENCODE's.
    """
    with open(path) as f:
        for line in f:
            if line.startswith("## HISTOGRAM"):
                break
        return np.loadtxt(f, skiprows=1, ndmin=2)


def nucleosomal_qc(data):
    """Compute ENCODE's fragment length checks.

    Parameters
    ----------
    data : numpy.ndarray
        Output of `read_histogram`.

    Returns
    -------
    list of str
        Lines of ENCODE's ``nucleosomal.qc``: metric, pass flag, and the
        value or the range a peak was missing from.
    """
    size, count = data[:, 0], data[:, 1]
    nfr = count[size < NFR_UPPER].sum()
    mono = count[(size > MONO_LOWER) & (size <= MONO_UPPER)].sum()
    lines = []
    for metric, value, lower in (
        ("Fraction of reads in NFR", nfr / count.sum(), 0.4),
        ("NFR / mono-nuc reads", nfr / mono, 2.5),
    ):
        ok = lower <= value
        lines.append(
            f"{metric}\t{ok}\t{value}\t"
            + ("OK" if ok else f"out of range [{lower}, inf]")
        )
    # find_peaks_cwt returns row indices, so ranges shift by the first size.
    peaks = find_peaks_cwt(count, np.array([25]))
    for metric, lower, upper in PEAK_RANGES:
        lower, upper = lower - size[0], upper - size[0]
        ok = any(lower <= p <= upper for p in peaks)
        lines.append(
            f"{metric}\t{ok}\t"
            + ("OK" if ok else f"Cannot find element in range [{lower}, {upper}]")
        )
    return lines


def plot_histogram(data, output_file):
    """Plot fragment counts per length, shading the NFR and mono-nucleosome ranges.

    Parameters
    ----------
    data : numpy.ndarray
        Output of `read_histogram`.
    output_file : str
        Plot SVG to write. The PNG goes beside it.
    """
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(data[:, 0], data[:, 1], width=1)
    colors = pk.PALETTES["plain"].cycle
    ax.axvspan(0, NFR_UPPER, color=colors[1], alpha=0.15, label="NFR")
    ax.axvspan(MONO_LOWER, MONO_UPPER, color=colors[2], alpha=0.15, label="Mono-nuc")
    ax.set(
        xlabel="Fragment length (bp)",
        ylabel="Fragments",
        title="Fragment length distribution",
        xlim=(0, 1000),
    )
    ax.legend(loc="upper right")
    pk.despine(ax)
    pk.save_figure(fig, output_file)
    plt.close(fig)


def main(argv=None):
    """Write ENCODE's fragment length checks and the distribution plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    data = read_histogram(args.histogram)
    with open(args.qc, "w") as f:
        f.writelines(line + "\n" for line in nucleosomal_qc(data))
    plot_histogram(data, args.plot)


if __name__ == "__main__":
    main()
