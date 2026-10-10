"""Compute library complexity (NRF, PBC1, PBC2) from BEDPE read pairs, as ENCODE does."""

import argparse
import json
import sys

import polars as pl

# bedtools bamtobed -bedpe columns ENCODE keys a fragment on.
KEY = {
    "column_1": "chrom1",
    "column_2": "start1",
    "column_4": "chrom2",
    "column_6": "end2",
    "column_9": "strand1",
    "column_10": "strand2",
}


def build_parser():
    """Build the pbc_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the pbc_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Compute NRF, PBC1 and PBC2.")
    p.add_argument("--mito", default="chrM", help="Mitochondrial contig to skip.")
    p.add_argument("--output", required=True, help="JSON to write.")
    return p


def library_complexity(bedpe, mito):
    """Count fragments by position and derive ENCODE's complexity metrics.

    A fragment is keyed on its first end's contig and start, its second end's
    contig and end, and both strands. Pairs whose first end is on ``mito``
    are skipped.

    Parameters
    ----------
    bedpe : polars.DataFrame
        Columns ``chrom1``, ``start1``, ``chrom2``, ``end2``, ``strand1`` and
        ``strand2``, one row per read pair.
    mito : str
        Mitochondrial contig.

    Returns
    -------
    dict of str to float
        Read pair totals and NRF, PBC1 and PBC2. PBC2 is -1 when no fragment
        is seen exactly twice, and a ratio is 0 when its denominator is 0.
    """
    counts = (
        bedpe.filter(pl.col("chrom1") != mito).group_by(list(KEY.values())).len()["len"]
    )
    total = int(counts.sum())
    distinct = counts.len()
    one = int((counts == 1).sum())
    two = int((counts == 2).sum())
    return {
        "total_read_pairs": total,
        "distinct_read_pairs": distinct,
        "one_read_pair": one,
        "two_read_pairs": two,
        "NRF": distinct / total if total else 0.0,
        "PBC1": one / distinct if distinct else 0.0,
        "PBC2": one / two if two else -1.0,
    }


def main(argv=None):
    """Read BEDPE from stdin and write the complexity metrics as JSON.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    bedpe = pl.read_csv(
        sys.stdin.buffer,
        separator="\t",
        has_header=False,
        columns=list(KEY),
        schema_overrides={"column_1": pl.String, "column_4": pl.String},
    ).rename(KEY)
    with open(args.output, "w") as f:
        json.dump(library_complexity(bedpe, args.mito), f, indent=4)


if __name__ == "__main__":
    main()
