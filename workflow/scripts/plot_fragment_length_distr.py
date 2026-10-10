"""Plot the fragment length distribution of a BAM and write nucleosome-range counts."""

import argparse
from collections import Counter

import matplotlib.pyplot as plt
import poikilos as pk
import pysam


def build_parser():
    """Build the plot_fragment_length_distr.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_fragment_length_distr.py flags.
    """
    p = argparse.ArgumentParser(
        description="Plot fragment length distribution and save statistics."
    )
    p.add_argument("--bam", required=True, help="Filtered BAM.")
    p.add_argument(
        "--plot", required=True, help="Distribution SVG to write (PNG beside it)."
    )
    p.add_argument("--output", required=True, help="Fragment statistics TXT to write.")
    return p


def calculate_fragment_length_distribution(bam_file):
    """Count the non-negative template lengths in a BAM.

    Parameters
    ----------
    bam_file : str
        BAM to read.

    Returns
    -------
    dict of int to int
        Count of reads per template length.
    """
    try:
        all_fragment_lengths = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            fragment_lengths = bam.fetch(until_eof=True)
            for fragment in fragment_lengths:
                if fragment.template_length < 0:
                    continue
                all_fragment_lengths.append(fragment.template_length)

        distribution = Counter(all_fragment_lengths)
        return dict(distribution)

    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def nucleosome_ranges(max_length):
    """Return the fragment length range of each nucleosome class.

    Parameters
    ----------
    max_length : int
        Upper end of the tri-nucleosomal range.

    Returns
    -------
    dict of str to tuple of int
        Inclusive (start, end) lengths per class label.
    """
    return {
        "Nucleosome-free (< 147bp)": (0, 147),
        "Mono-nucleosomal (147-294bp)": (147, 294),
        "Di-nucleosomal (294-441bp)": (294, 441),
        "Tri-nucleosomal (> 441bp)": (441, max_length),
    }


def range_count(distribution, start, end, max_length):
    """Count the fragments with a length from ``start`` to ``end``, inclusive.

    Parameters
    ----------
    distribution : dict of int to int
        Count of fragments per length.
    start : int
        Shortest length counted.
    end : int
        Longest length counted, capped at ``max_length``.
    max_length : int
        Longest length counted in any range.

    Returns
    -------
    int
        Fragments in the range.
    """
    return sum(
        distribution.get(length, 0)
        for length in range(start, min(end + 1, max_length + 1))
    )


def plot_fragment_length_distribution(distribution, output_file, max_length=1000):
    """Plot fragment counts per length, with the nucleosome ranges in the legend.

    Parameters
    ----------
    distribution : dict of int to int
        Count of fragments per length.
    output_file : str
        Plot SVG to write. The PNG goes beside it.
    max_length : int, default 1000
        Longest fragment length shown.
    """
    lengths = sorted(distribution.keys())
    counts = [distribution[length] for length in lengths]

    limited_lengths = [length for length in lengths if length <= max_length]
    limited_counts = [distribution[length] for length in limited_lengths]

    total_fragments = sum(counts)

    ranges = nucleosome_ranges(max_length)
    percentages = {
        label: (range_count(distribution, start, end, max_length) / total_fragments)
        * 100
        for label, (start, end) in ranges.items()
    }

    pk.use("plain")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(limited_lengths, limited_counts, width=1)
    ax.set(
        xlabel="Fragment Length (bp)",
        ylabel="Frequency",
        title="Fragment Length Distribution",
        xlim=(0, max_length),
    )

    for _, end in ranges.values():
        ax.axvline(x=end, color=pk.PALETTES["plain"].cycle[3], linestyle="--")

    ax.legend(
        [f"{label}: {percentages[label]:.1f}%" for label in percentages],
        loc="upper right",
    )
    pk.despine(ax)
    pk.save_figure(fig, output_file)
    plt.close(fig)
    print(f"Plot saved to {output_file}")


def save_fragment_statistics(distribution, output_file, max_length=1000):
    """Write the total fragment count and the count and share of each range.

    Parameters
    ----------
    distribution : dict of int to int
        Count of fragments per length.
    output_file : str
        Statistics TXT to write.
    max_length : int, default 1000
        Longest fragment length counted in a range.
    """
    total_fragments = sum(distribution.values())

    with open(output_file, "w") as f:
        f.write(f"Total fragments: {total_fragments}\n")
        for label, (start, end) in nucleosome_ranges(max_length).items():
            count = range_count(distribution, start, end, max_length)
            percentage = (count / total_fragments) * 100
            f.write(f"{label}: {count} fragments ({percentage:.2f}%)\n")
        print(f"Fragment statistics saved to {output_file}")


def main(argv=None):
    """Count fragment lengths, then write the plot and the statistics.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    try:
        distribution = calculate_fragment_length_distribution(args.bam)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise

    plot_fragment_length_distribution(distribution, args.plot)

    save_fragment_statistics(distribution, args.output)


if __name__ == "__main__":
    main()
