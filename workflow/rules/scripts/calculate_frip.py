import subprocess
import argparse

def calculate_fraction_reads_in_peaks(fragments_bed, peaks_bed):
    """
    Calculate the fraction of reads (fragments) that overlap peaks using bedtools.

    :param fragments_bed: Path to the BED file of fragments.
    :param peaks_bed: Path to the BED file of peaks.
    :return: Fraction of reads in peaks (float).
    """
    try:
        # Use bedtools intersect to find overlapping fragments and peaks
        cmd_overlap = f"bedtools intersect -u -a {fragments_bed} -b {peaks_bed} | wc -l"
        cmd_total = f"zcat {fragments_bed} | wc -l"
        
        # Count the number of overlapping fragments
        overlapping_count = int(subprocess.check_output(cmd_overlap, shell=True).strip())
        
        # Count the total number of fragments
        total_count = int(subprocess.check_output(cmd_total, shell=True).strip())
        
        # Calculate the fraction of reads in peaks
        fraction_in_peaks = overlapping_count / total_count if total_count > 0 else 0
        
        return fraction_in_peaks

    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools: {e.stderr}")
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(description="Calculate the fraction of reads in peaks.")
    parser.add_argument("fragments_bed", help="Path to the BED file of fragments")
    parser.add_argument("peaks_bed", help="Path to the BED file of peaks")
    parser.add_argument("output_file", help="Path to the output text file to write the fraction of reads in peaks")

    args = parser.parse_args()

    # Calculate the fraction of reads in peaks
    fraction = calculate_fraction_reads_in_peaks(args.fragments_bed, args.peaks_bed)

    # Write the result to the output file
    with open(args.output_file, "w") as f:
        f.write(f"Fraction of reads in peaks: {fraction}\n")
    
    print(f"Fraction of reads in peaks written to {args.output_file}")


if __name__ == "__main__":
    main()
