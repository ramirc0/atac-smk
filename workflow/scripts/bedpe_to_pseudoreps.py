"""Split a BEDPE into two pseudoreplicate tagAlign files by shuffling its rows."""

import argparse

import polars as pl

RANDOM_SEED = 61


def build_parser():
    """Build the bedpe_to_pseudoreps.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the bedpe_to_pseudoreps.py flags.
    """
    p = argparse.ArgumentParser(
        description="Split a BEDPE file into two pseudoreplicated BEDPE files."
    )
    p.add_argument("--bedpe", required=True, help="BEDPE to split.")
    p.add_argument("--pseudorep1", required=True, help="First tagAlign to write.")
    p.add_argument("--pseudorep2", required=True, help="Second tagAlign to write.")
    return p


def main(argv=None):
    """Shuffle the BEDPE rows, halve them and write one tagAlign per half.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    positive_columns = ["column_1", "column_2", "column_3", "column_9"]
    negative_columns = ["column_4", "column_5", "column_6", "column_10"]

    print("Reading BEDPE file...")
    bedpe_df = pl.read_csv(
        args.bedpe, separator="\t", has_header=False, comment_prefix="#"
    )

    print("Shuffling BEDPE file...")
    shuffled_df = bedpe_df.sample(fraction=1.0, shuffle=True, seed=RANDOM_SEED)
    length = bedpe_df.shape[0]

    print("Splitting BEDPE file into pseudoreplicated BEDPE files...")
    pseudorep1_df = shuffled_df.slice(0, length // 2)
    pseudorep2_df = shuffled_df.slice(length // 2, length)

    print("Converting to tagalign format...")
    pseudorep1_df_pos = pseudorep1_df.drop(negative_columns)
    pseudorep1_df_neg = pseudorep2_df.drop(positive_columns)
    default_names = [
        f"column_{i}" for i in range(1, len(pseudorep1_df_neg.columns) + 1)
    ]
    pseudorep1_df_pos.columns = default_names
    pseudorep1_df_neg.columns = default_names
    pseudorep1_df_concat = pl.concat([pseudorep1_df_pos, pseudorep1_df_neg])

    pseudorep2_df_pos = pseudorep2_df.drop(positive_columns)
    pseudorep2_df_neg = pseudorep2_df.drop(negative_columns)
    pseudorep2_df_pos.columns = default_names
    pseudorep2_df_neg.columns = default_names
    pseudorep2_df_concat = pl.concat([pseudorep2_df_pos, pseudorep2_df_neg])

    print("Writing pseudoreplicated tagalign files...")
    pseudorep1_df_concat.write_csv(
        args.pseudorep1, include_header=False, separator="\t"
    )
    pseudorep2_df_concat.write_csv(
        args.pseudorep2, include_header=False, separator="\t"
    )

    print("Done!")


if __name__ == "__main__":
    main()
