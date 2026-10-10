import subprocess
from collections import Counter
import matplotlib.pyplot as plt
import argparse
import pysam


def calculate_fragment_length_distribution(bam_file):
    """ Calculate the fragment length distribution from a BAM file using samtools.

    :param bam_file: Path to the input BAM file.
    :param threads: Number of threads for samtools (default: 1).
    :return: A dictionary where keys are fragment lengths and values are their frequencies.
    """

    try:
        # Step 1: Extract fragment lengths using pysam
        all_fragment_lengths = []
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            fragment_lengths = bam.fetch(until_eof=True)
            for fragment in fragment_lengths:
                if fragment.template_length < 0:
                    continue
                all_fragment_lengths.append(fragment.template_length)
            
        distribution = Counter(all_fragment_lengths)
        return dict(distribution)

    except subprocess.CalledProcessError as e:
        print(f"Error running samtools or awk: {e.stderr}")
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def plot_fragment_length_distribution(distribution, output_file, max_length=1000):
    """
    Create a bar plot of fragment length distribution with nucleosome-associated regions,
    including the percentage of fragments in each range in the legend.
    
    :param distribution: Dictionary where keys are fragment lengths and values are counts.
    :param output_file: Path to save the plot as an image file.
    :param max_length: Maximum fragment length to display on the x-axis (default: 1000 bp).
    """
    # Sort the fragment lengths
    lengths = sorted(distribution.keys())
    counts = [distribution[length] for length in lengths]

    # Limit the range for visualization
    limited_lengths = [l for l in lengths if l <= max_length]
    limited_counts = [distribution[l] for l in limited_lengths]

    # Calculate total number of fragments
    total_fragments = sum(counts)

    # Define nucleosome regions and calculate percentages
    nucleosome_ranges = {
        "Nucleosome-free (< 147bp)": (0, 147),
        "Mono-nucleosomal (147-294bp)": (147, 294),
        "Di-nucleosomal (294-441bp)": (294, 441),
        "Tri-nucleosomal (> 441bp)": (441, max_length),
    }

    percentages = {}
    for label, (start, end) in nucleosome_ranges.items():
        # Sum counts within the range
        count = sum(distribution.get(length, 0) for length in range(start, min(end + 1, max_length + 1)))
        percentages[label] = (count / total_fragments) * 100

    # Create the bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(limited_lengths, limited_counts, width=1, color="skyblue", edgecolor="black", alpha=0.7)
    plt.xlabel("Fragment Length (bp)", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title("Fragment Length Distribution", fontsize=16, pad=15)
    plt.xlim(0, max_length)
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Add vertical lines and legend
    for label, (start, end) in nucleosome_ranges.items():
        plt.axvline(x=end, color="red", linestyle="--", linewidth=1.5)

    plt.legend(
        [f"{label}: {percentages[label]:.1f}%" for label in percentages],
        loc="upper right",
        fontsize=10,
    )

    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {output_file}")


def save_fragment_statistics(distribution, output_file, max_length=1000):
    """
    Save the number of fragments and percentages for each range to a text file.

    :param distribution: Dictionary where keys are fragment lengths and values are counts.
    :param output_file: Path to save the statistics as a text file.
    :param max_length: Maximum fragment length considered (default: 1000 bp).
    """
    total_fragments = sum(distribution.values())

    # TODO: Make this JSON smh

    # Define nucleosome regions
    nucleosome_ranges = {
        "Nucleosome-free (< 147bp)": (0, 147),
        "Mono-nucleosomal (147-294bp)": (147, 294),
        "Di-nucleosomal (294-441bp)": (294, 441),
        "Tri-nucleosomal (> 441bp)": (441, max_length),
    }

    with open(output_file, "w") as f:
        f.write(f"Total fragments: {total_fragments}\n")
        for label, (start, end) in nucleosome_ranges.items():
            count = sum(distribution.get(length, 0) for length in range(start, min(end + 1, max_length + 1)))
            percentage = (count / total_fragments) * 100
            f.write(f"{label}: {count} fragments ({percentage:.2f}%)\n")
        print(f"Fragment statistics saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Plot fragment length distribution and save statistics.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("threads", type=int, help="Number of threads for samtools")
    parser.add_argument("output_png", help="Output PNG file for the plot")
    parser.add_argument("output_txt", help="Output TXT file for fragment statistics")
    args = parser.parse_args()

    # Step 1: Calculate fragment length distribution
    try:
        distribution = calculate_fragment_length_distribution(args.bam_file)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise

    # Step 2: Plot fragment length distribution
    plot_fragment_length_distribution(distribution, args.output_png)

    # Step 3: Save fragment statistics
    save_fragment_statistics(distribution, args.output_txt)


if __name__ == "__main__":
    main()
