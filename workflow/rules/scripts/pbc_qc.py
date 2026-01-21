import oxbow
import polars as pl
import json
import argparse

def main(args):
    bed_cols = ["chrom", "start", "end", "name", "score", "strand"]
    agg_cols = ["chrom", "start", "end", "strand"]
    
    bam_file = oxbow.from_bam(args.bam).to_polars(lazy=True)
    bed_lf = bam_file.with_columns(
        pl.col("rname").alias("chrom"),
        pl.min_horizontal("pos", "pnext").alias("start"),
        (pl.max_horizontal("pos", "pnext") + pl.col("tlen")).alias("end"),
        pl.col("mapq").alias("score"),
        pl.col("qname").alias("name"),
        pl.when(pl.col("flag") & 16 == 0)
        .then(pl.lit("+"))
        .otherwise(pl.lit("-"))
        .alias("strand")
    ).select(bed_cols + ["flag"]).filter(
        ~((pl.col("chrom") == "chrM") | (pl.col("chrom").cat.starts_with("chrUn"))), # Filter mitochondrial and reads mapped to contigs/scaffolds
        ((pl.col("flag") & 1804 == 0) & (pl.col("flag") & 2 != 0)), # filtering for (same as filter.smk)
        (pl.col("score") >= 30)
    ).group_by("name").agg(
        [pl.col(col).filter(pl.col("flag") & 64 == 0).first().alias(col) for col in agg_cols] + \
        [pl.col(col).filter(pl.col("flag") & 128 == 0).first().alias(col + "_") for col in agg_cols] + \
        [pl.col("score").first()]
    ).select(
        ["chrom", "start", "end", "chrom_", "start_", "end_", "name", "score", "strand", "strand_"]
    ).sort("name")

    bed_df = bed_lf.collect()
    uniq_counts = bed_df.with_columns(
        pl.concat_str(["chrom", "start", "end", "chrom_", "start_", "end_"]).alias("uniq")
    )["uniq"].value_counts()

    no_total_reads = uniq_counts["count"].sum()
    no_distict_reads = uniq_counts.height
    no_unique_reads = uniq_counts.filter(pl.col("count") == 1).height
    no_reads_seen_twice = uniq_counts.filter(pl.col("count") == 2).height

    print(f"Total reads: {no_total_reads}")
    print(f"Distinct reads: {no_distict_reads}")
    print(f"Unique reads: {no_unique_reads}")
    print(f"Reads seen twice: {no_reads_seen_twice}")

    PBC1 = no_total_reads / no_distict_reads
    if no_reads_seen_twice != 0:
        PBC2 = no_total_reads / no_reads_seen_twice
    else:
        PBC2 = float('inf')
    NRF = no_distict_reads / no_total_reads

    print(f"PBC1: {PBC1}")
    print(f"PBC2: {PBC2}")
    print(f"NRF: {NRF}")

    with open(str(args.pbc_json), "w") as f:
        json.dump(dict(PBC1=PBC1, PBC2=PBC2, NRF=NRF), f, indent=4)

    print("Wrote QC metrics to ", str(args.pbc_json))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("pbc_json")

    args = parser.parse_args()

    main(args)
    