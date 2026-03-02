import polars as pl

def main():
    ( 
        pl
        .read_csv("./MOHD.Phase_0_Samples.tsv", separator="\t")
        .filter(pl.col("data_modality") == "ATAC-seq")
        .with_columns(
            pl.col("entity:dataset_id").alias("mohd_accession"),
        )
        .select(
            "mohd_accession",
            pl.format("/data/projects/mohd/data/Molecular/2_ATAC/{}/{}_R1.fastq.gz", "mohd_accession", "mohd_accession").alias("r1"),
            pl.format("/data/projects/mohd/data/Molecular/2_ATAC/{}/{}_R2.fastq.gz", "mohd_accession", "mohd_accession").alias("r2"),
        )
        .with_columns(
            pl.col("mohd_accession").str.strip_chars_start("MOHD_")
        )
        .sort("mohd_accession")
        .write_csv("./MOHD.sample_map.tsv", separator="\t")
    )

if __name__ == "__main__":
    main()
