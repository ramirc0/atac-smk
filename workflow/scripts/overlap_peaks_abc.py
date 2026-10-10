import pandas as pd
import bioframe as bf
import argparse
from pathlib import Path

def filter_50percent_overlap(df1, df2):
    """Find rows in df1 that overlap with df2 by at least 50%"""
    if df1.empty or df2.empty:
        return df1.iloc[:0]
    
    overlaps = bf.overlap(df1, df2, return_overlap=True, how='inner')
    
    if overlaps.empty:
        return df1.iloc[:0]
    
    # Calculate overlap fractions
    overlaps['len1'] = overlaps['end'] - overlaps['start']
    overlaps['len2'] = overlaps['end_'] - overlaps['start_']
    overlaps['overlap_len'] = overlaps['overlap_end'] - overlaps['overlap_start']
    overlaps['frac1'] = overlaps['overlap_len'] / overlaps['len1']
    overlaps['frac2'] = overlaps['overlap_len'] / overlaps['len2']
    
    # Filter for 50% overlap
    mask = (overlaps['frac1'] >= 0.5) | (overlaps['frac2'] >= 0.5)
    filtered = overlaps[mask]
    
    # Return the original rows from df1 that passed the filter
    passing_rows = filtered[['chrom', 'start', 'end', 'name', 'score', 'strand', 
                             'signalValue', 'pValue', 'qValue', 'peak']].drop_duplicates()
    
    return passing_rows

def main(args):
    all_fragments = args.all_fragments
    subset1 = args.subset1
    subset2 = args.subset2
    output = args.output
    
    peak_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                 'signalValue', 'pValue', 'qValue', 'peak']
    
    all_fragments = pd.read_csv(all_fragments, sep='\t', names=peak_cols, header=None)
    subset1 = pd.read_csv(subset1, sep='\t', names=peak_cols, header=None)
    subset2 = pd.read_csv(subset2, sep='\t', names=peak_cols, header=None)
    
    # Find peaks that overlap subset1 by ≥50%
    passed_peaks = filter_50percent_overlap(all_fragments, subset1)
    
    if passed_peaks.empty:
        raise ValueError("No peaks that overlap subset 1 by ≥50%")
    
    # Find peaks that also overlap subset2 by ≥50%
    final_peaks = filter_50percent_overlap(passed_peaks, subset2)
    
    if final_peaks.empty:
        raise ValueError("No peaks that overlap both subsets by ≥50%")
    
    final_peaks = final_peaks.drop_duplicates().sort_values(['chrom', 'start'])
    final_peaks.to_csv(output, sep='\t', header=False, index=False)
    print(f"Found {len(final_peaks)} peaks that overlap both subsets by ≥50%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("all_fragments", type=Path, help="Path to all fragments narrowPeak file")
    parser.add_argument("subset1", type=Path, help="Path to subset 1 narrowPeak file")
    parser.add_argument("subset2", type=Path, help="Path to subset 2 narrowPeak file")
    parser.add_argument("output", type=Path, help="Path to output overlap peaks file")
    args = parser.parse_args()
    main(args)