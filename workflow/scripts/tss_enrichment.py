import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt
import json
import os

def ensure_bam_index(bam_file):
    """Ensure the BAM file has an index. If the index is missing or outdated, create one."""
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
    """Extract signal from a BAM file centered on the cut site within a ±window around TSSs."""
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
    """Calculate the TSS enrichment score from the signal matrix."""
    X_norm = matrix / np.nanmean(np.concatenate([matrix[:, :100], matrix[:, -100:]], axis=1)) # This would represent the background signal presumably?
    enrichment = np.max(np.nanmean(X_norm, axis=0))
    return enrichment

def plot_tss_enrichment(matrix, enrichment_score, output_plot):
    """Plot TSS enrichment profile and save the figure."""
    mean_signal = np.nanmean(matrix, axis=0)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(np.arange(len(mean_signal)), mean_signal, c="red")
    ax.axvline(x=2000, color="black", linestyle="--", label="TSS")
    ax.set_xlim([0, 4000])
    ax.set_xticks([0, 2000, 4000])
    ax.set_xticklabels(["-2kb", "TSS", "+2kb"])
    ax.set_xlabel("Position relative to TSS", fontsize=12)
    ax.set_ylabel("Average Signal", fontsize=12)
    ax.set_title(f"TSS Enrichment (Score: {enrichment_score:.2f})", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    print(f"TSS enrichment plot saved to {output_plot}")

def main():
    parser = argparse.ArgumentParser(description="Calculate and plot TSS enrichment.")
    parser.add_argument("bam_file", help="Path to the BAM file.")
    parser.add_argument("tss_bed", help="Path to the TSS BED file.")
    parser.add_argument("output_png", help="Path to save the TSS enrichment plot (PNG).")
    parser.add_argument("output_txt", help="Path to save the TSS enrichment score (TXT).")
    parser.add_argument("--save-npy", help="Optional: Path to save the signal matrix as a .npy file.")
    args = parser.parse_args()
    
    ensure_bam_index(args.bam_file)
    
    tss_df = pd.read_csv(
        args.tss_bed, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"]
    )
    
    signal_matrix = extract_signal_around_tss_with_cut_site(args.bam_file, tss_df)
    enrichment_score = calculate_tss_enrichment(signal_matrix)
    
    with open(args.output_txt, "w") as f:
        json.dump(dict(tss_enrichment_score=enrichment_score), f, indent=4)
    print(f"TSS enrichment score saved to {args.output_txt}")
    
    plot_tss_enrichment(signal_matrix, enrichment_score, args.output_png)
    
    if args.save_npy:
        np.save(args.save_npy, signal_matrix)
        print(f"Signal matrix saved to {args.save_npy}")

if __name__ == "__main__":
    main()
