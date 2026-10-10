"""Compute the fraction of tagAlign reads that fall in peaks (FRiP)."""

import argparse
import subprocess


def build_parser():
    """Build the calculate_frip.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the calculate_frip.py flags.
    """
    p = argparse.ArgumentParser(description="Calculate the fraction of reads in peaks.")
    p.add_argument("--tagalign", required=True, help="Gzipped tagAlign of reads.")
    p.add_argument("--peaks", required=True, help="Peaks BED.")
    p.add_argument("--output", required=True, help="Text file for the FRiP line.")
    return p


def calculate_fraction_reads_in_peaks(fragments_bed, peaks_bed):
    """Compute the fraction of reads that overlap a peak, using bedtools.

    Parameters
    ----------
    fragments_bed : str
        Gzipped BED of reads.
    peaks_bed : str
        BED of peaks.

    Returns
    -------
    float
        Reads overlapping a peak over all reads, or 0 if there are no reads.
    """
    try:
        cmd_overlap = f"bedtools intersect -u -a {fragments_bed} -b {peaks_bed} | wc -l"
        cmd_total = f"zcat {fragments_bed} | wc -l"

        overlapping_count = int(
            subprocess.check_output(cmd_overlap, shell=True).strip()
        )
        total_count = int(subprocess.check_output(cmd_total, shell=True).strip())

        fraction_in_peaks = overlapping_count / total_count if total_count > 0 else 0

        return fraction_in_peaks

    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools: {e.stderr}")
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def main(argv=None):
    """Compute FRiP and write it as one line of text.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    fraction = calculate_fraction_reads_in_peaks(args.tagalign, args.peaks)

    with open(args.output, "w") as f:
        f.write(f"Fraction of reads in peaks: {fraction}\n")

    print(f"Fraction of reads in peaks written to {args.output}")


if __name__ == "__main__":
    main()
