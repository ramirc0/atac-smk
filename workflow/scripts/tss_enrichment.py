"""Score and plot TSS enrichment with ENCODE's metaseq algorithm.

Port of ENCODE's ``encode_task_tss_enrich.py`` and the metaseq 0.5.6 code it
calls, quirks included, so scores compare with ENCODE's. Every aligned base
of every read counts, after the read moves half a read length toward its 5'
end. Reads are fetched from the window shrunk by that half length, so near
the window edges some reads are missed. metaseq's CIGAR walk also advances
past insertions and soft clips as if they were on the reference.
"""

import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk
import polars as pl
import pysam
from zlib_ng import gzip_ng

BINS = 400
EDGE = 2000
EDGE_BINS = int(100 / (2 * EDGE / BINS))


def build_parser():
    """Build the tss_enrichment.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the tss_enrichment.py flags.
    """
    p = argparse.ArgumentParser(description="Score and plot TSS enrichment.")
    p.add_argument("--bam", required=True, help="Filtered, indexed BAM.")
    p.add_argument("--tss", required=True, help="TSS BED, gzipped or not.")
    p.add_argument("--chromsizes", required=True, help="Contig sizes.")
    p.add_argument("--fastq", required=True, help="Trimmed R1 FASTQ, for read length.")
    p.add_argument("--qc", required=True, help="tss_enrich.qc to write.")
    p.add_argument(
        "--plot", required=True, help="Profile SVG to write (PNG beside it)."
    )
    return p


def read_length(fastq):
    """Return the longest of the first million reads, ENCODE's read length.

    Parameters
    ----------
    fastq : str
        Gzipped FASTQ.

    Returns
    -------
    int
        Longest sequence among the first million reads.
    """
    longest = 0
    with gzip_ng.open(fastq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                longest = max(longest, len(line.strip()))
            if i == 4_000_000:
                break
    return longest


def coverage(bam, chrom, start, stop, shift):
    """Return metaseq's per-base coverage of one window.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed BAM.
    chrom : str
        Window contig.
    start, stop : int
        Window bounds, 0-based, half-open.
    shift : int
        Negative ``shift_width`` of metaseq: minus half the read length.

    Returns
    -------
    numpy.ndarray
        Count of shifted aligned bases at each window position.
    """
    size = stop - start
    steps = np.zeros(size + 1)
    for read in bam.fetch(chrom, max(start - shift, 0), stop + shift):
        move = -shift if read.is_reverse else shift
        end = read.reference_start
        for op, length in read.cigartuples:
            block, end = end, end + length
            if op != 0:
                continue
            i = max(block + move - start, 0)
            j = min(end + move - start, size)
            if i >= size or j < 0:
                continue
            steps[i] += 1
            steps[j] -= 1
    return np.cumsum(steps[:-1])


def tss_matrix(bam_file, tss, chromsizes, read_len):
    """Build metaseq's TSS-by-bin coverage matrix.

    Parameters
    ----------
    bam_file : str
        Indexed BAM.
    tss : polars.DataFrame
        ``chrom``, ``start``, ``end`` and ``strand`` per TSS.
    chromsizes : dict of str to int
        Contig sizes, to clip the windows.
    read_len : int
        ENCODE's read length.

    Returns
    -------
    numpy.ndarray
        One row of ``BINS`` point samples per TSS, 5' to 3' on its strand.
    """
    # Python 2 floor division, as ENCODE ran it: -read_len/2.
    shift = math.floor(-read_len / 2)
    rows = []
    with pysam.AlignmentFile(bam_file) as bam:
        for chrom, start, end, strand in tss.iter_rows():
            start = max(start - EDGE, 0)
            stop = min(end + EDGE, chromsizes[chrom])
            profile = np.interp(
                np.linspace(start, stop - 1, BINS),
                np.arange(start, stop),
                coverage(bam, chrom, start, stop, shift),
            )
            rows.append(profile[::-1] if strand == "-" else profile)
    return np.array(rows)


def enrichment(matrix):
    """Normalize to the flanks and score, as ENCODE's Greenleaf normalization.

    Parameters
    ----------
    matrix : numpy.ndarray
        Output of `tss_matrix`.

    Returns
    -------
    numpy.ndarray
        Mean normalized profile.
    numpy.float64
        TSS enrichment: the profile's maximum.
    """
    means = matrix.mean(axis=0)
    noise = (sum(means[:EDGE_BINS]) + sum(means[-EDGE_BINS:])) / (2 * EDGE_BINS)
    profile = (matrix / noise).mean(axis=0)
    return profile, max(profile)


def plot_profile(profile, score, output_file):
    """Plot the normalized mean profile around the TSS.

    Parameters
    ----------
    profile : numpy.ndarray
        Mean normalized profile.
    score : float
        TSS enrichment, shown in the title.
    output_file : str
        Plot SVG to write. The PNG goes beside it.
    """
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(np.linspace(-EDGE, EDGE, BINS), profile)
    ax.axvline(0, color=pk.PALETTES["plain"].foreground, linestyle=":")
    ax.set(
        xlabel="Distance from TSS (bp)",
        ylabel="TSS enrichment",
        title=f"TSS enrichment {score:.2f}",
    )
    pk.despine(ax)
    pk.save_figure(fig, output_file)
    plt.close(fig)


def main(argv=None):
    """Write the TSS enrichment score and the profile plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    tss = pl.read_csv(
        args.tss,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 5],
        new_columns=["chrom", "start", "end", "strand"],
    )
    chromsizes = dict(
        pl.read_csv(args.chromsizes, separator="\t", has_header=False).iter_rows()
    )
    matrix = tss_matrix(args.bam, tss, chromsizes, read_length(args.fastq))
    profile, score = enrichment(matrix)
    with open(args.qc, "w") as f:
        f.write(f"{score}")
    plot_profile(profile, score, args.plot)


if __name__ == "__main__":
    main()
