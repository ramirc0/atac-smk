"""Score and plot Tn5 cut-site enrichment in a ±2 kb window around TSSs."""

import argparse
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import poikilos as pk
import pysam
from tqdm import tqdm


def build_parser():
    """Build the tss_enrichment.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tss_enrichment.py flags.
    """
    p = argparse.ArgumentParser(description="Calculate and plot TSS enrichment.")
    p.add_argument("--bam", required=True, help="Filtered BAM.")
    p.add_argument("--tss", required=True, help="TSS BED.")
    p.add_argument(
        "--plot", required=True, help="Enrichment SVG to write (PNG beside it)."
    )
    p.add_argument("--output", required=True, help="Score JSON to write.")
    p.add_argument("--save-npy", help="Optional .npy to write the signal matrix to.")
    return p


def ensure_bam_index(bam_file):
    """Index the BAM if its index is missing or older than the BAM.

    Parameters
    ----------
    bam_file : str
        BAM to index.
    """
    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        print(f"Index for {bam_file} not found. Creating index...")
        pysam.index(bam_file)
        print(f"Index created for {bam_file}.")
    else:
        bam_mtime = os.path.getmtime(bam_file)
        index_mtime = os.path.getmtime(index_file)
        if index_mtime < bam_mtime:
            print(f"Index for {bam_file} is older than the BAM file. Re-indexing...")
            pysam.index(bam_file)
            print(f"Index updated for {bam_file}.")
        else:
            print(f"Index for {bam_file} is present and up-to-date.")


def extract_signal_around_tss_with_cut_site(bam_file, tss_df, window=2000):
    """Count Tn5 cut sites per position in a ±``window`` bp window around each TSS.

    Parameters
    ----------
    bam_file : str
        Indexed BAM.
    tss_df : pandas.DataFrame
        TSSs with ``chrom``, ``start`` and ``strand`` columns.
    window : int, default 2000
        Bases on each side of the TSS.

    Returns
    -------
    numpy.ndarray
        One row per TSS, oriented 5' to 3' on the TSS strand.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    matrix = []

    for _, tss in tqdm(tss_df.iterrows(), total=len(tss_df), desc="Extracting signal"):
        chrom = tss["chrom"]
        start = max(0, tss["start"] - window)
        end = tss["start"] + window
        strand = tss["strand"]

        coverage = np.zeros(2 * window + 1)
        for pileup_column in bam.pileup(chrom, start, end, truncate=True):
            cut_site_pos = None
            for pileup_read in pileup_column.pileups:
                if pileup_read.alignment.is_reverse:
                    cut_site_pos = pileup_read.alignment.reference_end - 5
                else:
                    cut_site_pos = pileup_read.alignment.reference_start + 4
                if cut_site_pos is not None and start <= cut_site_pos < end:
                    coverage[cut_site_pos - start] += 1
        if strand == "-":
            coverage = coverage[::-1]
        matrix.append(coverage)

    bam.close()
    return np.array(matrix)


def calculate_tss_enrichment(matrix):
    """Compute the peak of the mean signal over the background-scaled TSS profiles.

    Parameters
    ----------
    matrix : numpy.ndarray
        Cut-site counts, one row per TSS.

    Returns
    -------
    numpy.float64
        TSS enrichment score.
    """
    # The outer 100 bp on each side of the window are the background.
    X_norm = matrix / np.nanmean(
        np.concatenate([matrix[:, :100], matrix[:, -100:]], axis=1)
    )
    enrichment = np.max(np.nanmean(X_norm, axis=0))
    return enrichment


def plot_tss_enrichment(matrix, enrichment_score, output_plot):
    """Plot the mean cut-site profile around the TSS.

    Parameters
    ----------
    matrix : numpy.ndarray
        Cut-site counts, one row per TSS.
    enrichment_score : float
        Score shown in the title.
    output_plot : str
        Plot SVG to write. The PNG goes beside it.
    """
    mean_signal = np.nanmean(matrix, axis=0)

    pk.use("plain")
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(np.arange(len(mean_signal)), mean_signal)
    ax.axvline(
        x=2000, color=pk.PALETTES["plain"].foreground, linestyle="--", label="TSS"
    )
    ax.set_xlim([0, 4000])
    ax.set_xticks([0, 2000, 4000], labels=["-2kb", "TSS", "+2kb"])
    ax.set(
        xlabel="Position relative to TSS",
        ylabel="Average Signal",
        title=f"TSS Enrichment (Score: {enrichment_score:.2f})",
    )
    ax.legend()
    pk.despine(ax)
    pk.save_figure(fig, output_plot)
    plt.close(fig)
    print(f"TSS enrichment plot saved to {output_plot}")


def main(argv=None):
    """Score TSS enrichment, then write the score JSON and the profile plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    ensure_bam_index(args.bam)

    tss_df = pd.read_csv(
        args.tss,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )

    signal_matrix = extract_signal_around_tss_with_cut_site(args.bam, tss_df)
    enrichment_score = calculate_tss_enrichment(signal_matrix)

    with open(args.output, "w") as f:
        json.dump({"tss_enrichment_score": enrichment_score}, f, indent=4)
    print(f"TSS enrichment score saved to {args.output}")

    plot_tss_enrichment(signal_matrix, enrichment_score, args.plot)

    if args.save_npy:
        np.save(args.save_npy, signal_matrix)
        print(f"Signal matrix saved to {args.save_npy}")


if __name__ == "__main__":
    main()
