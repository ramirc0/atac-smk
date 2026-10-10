#!/usr/bin/env python

import polars as pl
import argparse

RANDOM_SEED = 61

def main(args):
    positive_columns = ["column_1", "column_2", "column_3", "column_9"]
    negative_columns = ["column_4", "column_5", "column_6", "column_10"]
    bedpe_file = args.bedpe_file
    pseudorep1_file = args.pseudorep1_file
    pseudorep2_file = args.pseudorep2_file
    
    print("Reading BEDPE file...")
    bedpe_df = pl.read_csv(bedpe_file, separator="\t", has_header=False, comment_prefix="#")
    
    print("Shuffling BEDPE file...")
    shuffled_df = bedpe_df.sample(fraction=1.0, shuffle=True, seed=RANDOM_SEED) 
    length = bedpe_df.shape[0]
    
    print("Splitting BEDPE file into pseudoreplicated BEDPE files...")
    pseudorep1_df = shuffled_df.slice(0, length//2)
    pseudorep2_df = shuffled_df.slice(length//2, length)

    print("Converting to tagalign format...")
    pseudorep1_df_pos = pseudorep1_df.drop(negative_columns)
    pseudorep1_df_neg = pseudorep2_df.drop(positive_columns)
    default_names = [f"column_{i}" for i in range(1, len(pseudorep1_df_neg.columns) + 1)]
    pseudorep1_df_pos.columns = default_names
    pseudorep1_df_neg.columns = default_names
    pseudorep1_df_concat = pl.concat([pseudorep1_df_pos, pseudorep1_df_neg])

    pseudorep2_df_pos = pseudorep2_df.drop(positive_columns)
    pseudorep2_df_neg = pseudorep2_df.drop(negative_columns)
    pseudorep2_df_pos.columns = default_names
    pseudorep2_df_neg.columns = default_names
    pseudorep2_df_concat = pl.concat([pseudorep2_df_pos, pseudorep2_df_neg])

    print("Writing pseudoreplicated tagalign files...")
    pseudorep1_df_concat.write_csv(pseudorep1_file, include_header=False, separator="\t")
    pseudorep2_df_concat.write_csv(pseudorep2_file, include_header=False, separator="\t")

    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a BEDPE file into two pseudoreplicated BEDPE files.")
    parser.add_argument("bedpe_file", help="Input BEDPE file")
    parser.add_argument("pseudorep1_file", help="Output pseudoreplicated BEDPE file 1")
    parser.add_argument("pseudorep2_file", help="Output pseudoreplicated BEDPE file 2")
    args = parser.parse_args()
    main(args)
