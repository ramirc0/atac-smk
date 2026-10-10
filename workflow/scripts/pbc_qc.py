"""Compute library complexity (PBC1, PBC2, NRF) from an aligned BAM."""

import argparse
import json

import oxbow
import polars as pl


def build_parser():
    """Build the pbc_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the pbc_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Compute PBC1, PBC2 and NRF.")
    p.add_argument("--bam", required=True, help="Aligned BAM.")
    p.add_argument(
        "--output", required=True, help="JSON of PBC1, PBC2 and NRF to write."
    )
    return p


def main(argv=None):
    """Count read pairs by position and write the library complexity metrics.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    bed_cols = ["chrom", "start", "end", "name", "score", "strand"]
    agg_cols = ["chrom", "start", "end", "strand"]

    bam_file = oxbow.from_bam(args.bam).to_polars(lazy=True)
    assert isinstance(bam_file, pl.LazyFrame)
    bed_lf = (
        bam_file.with_columns(
            pl.col("rname").alias("chrom"),
            pl.col("pos").alias("start"),
            (pl.col("pos") + pl.col("seq").str.len_chars()).alias("end"),
            pl.col("mapq").alias("score"),
            pl.col("qname").alias("name"),
            pl.when((pl.col("flag") & 16) == 0)
            .then(pl.lit("+"))
            .otherwise(pl.lit("-"))
            .alias("strand"),
        )
        .select(bed_cols + ["flag"])
        .filter(
            ~((pl.col("chrom") == "chrM") | (pl.col("chrom").cat.starts_with("chrUn"))),
            # Hardcodes the filter rule's default samtools flags.
            (((pl.col("flag") & 1804) == 0) & ((pl.col("flag") & 2) != 0)),
            (pl.col("score") >= 30),
        )
        .group_by("name")
        .agg(
            [
                pl.col(col).filter((pl.col("flag") & 64) == 64).first().alias(col)
                for col in agg_cols
            ]
            + [
                pl.col(col)
                .filter((pl.col("flag") & 128) == 128)
                .first()
                .alias(col + "_")
                for col in agg_cols
            ]
            + [pl.col("score").first()]
        )
        .select(
            [
                "chrom",
                "start",
                "end",
                "chrom_",
                "start_",
                "end_",
                "name",
                "score",
                "strand",
                "strand_",
            ]
        )
        .sort("name")
    )

    bed_df = bed_lf.collect()
    uniq_counts = bed_df.with_columns(
        pl.concat_str(["chrom", "start", "end", "chrom_", "start_", "end_"]).alias(
            "uniq"
        )
    )["uniq"].value_counts()

    no_total_reads = uniq_counts["count"].sum()
    no_distict_reads = uniq_counts.height
    no_unique_reads = uniq_counts.filter(pl.col("count") == 1).height
    no_reads_seen_twice = uniq_counts.filter(pl.col("count") == 2).height

    print(f"Total reads: {no_total_reads}")
    print(f"Distinct reads: {no_distict_reads}")
    print(f"Unique reads: {no_unique_reads}")
    print(f"Reads seen twice: {no_reads_seen_twice}")

    PBC1 = no_unique_reads / no_distict_reads
    if no_reads_seen_twice != 0:
        PBC2 = no_unique_reads / no_reads_seen_twice
    else:
        PBC2 = float("inf")
    NRF = no_distict_reads / no_total_reads

    print(f"PBC1: {PBC1}")
    print(f"PBC2: {PBC2}")
    print(f"NRF: {NRF}")

    with open(args.output, "w") as f:
        json.dump({"PBC1": PBC1, "PBC2": PBC2, "NRF": NRF}, f, indent=4)

    print("Wrote QC metrics to ", args.output)


if __name__ == "__main__":
    main()
