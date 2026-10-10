# ATAC-seq Snakemake Pipeline

## Overview

Here we present a [Snakemake](https://github.com/snakemake/snakemake.git) pipeline for ATAC-seq data processing, peak calling, quality control, and reproducibility analysis across multiple samples.

This pipeline implements the [ENCODE ATAC-seq processing standards](https://www.encodeproject.org/atac-seq/), covering adapter trimming, alignment, filtering, signal track generation, peak calling, and irreproducible discovery rate (IDR) analysis.

## Example Usage

### Prerequisites

To run this pipeline, you will need paired-end FASTQ files from ATAC-seq experiments. A sample map TSV file must be prepared that maps sample accessions to their corresponding FASTQ paths (see [Example Configfile](#example-configfile)).

You will also need the following reference files placed in a `resources/` directory:

- GRCh38 reference genome FASTA (no-alt analysis set)
- Reference genome index (`.fai`) and chromosome sizes
- TSS reference BED file (e.g., `ENCFF493CCB.bed`)

### Dependencies

This pipeline requires the following dependencies:

- [`pixi`](https://pixi.sh)
- (Optional) [`slurm`](https://slurm.schedmd.com/quickstart.html)

To start, set up the launcher environment:

```bash
pixi shell
```

> [!NOTE]
> Rules run in two pixi environments from `workflow/envs/pixi.toml`: `default` (for most rules) and `idr`, since IDR requires a pinned Python 3.9 / NumPy 1.19 stack.

### Example Configfile

For a complete example, see [config/config.yml](config/config.yml).

#### Input

```yaml
input:
  sample_map: "MOHD.sample_map.tsv"
  tmp_dir: ".tmp"
```

- `input.sample_map`: Path to your sample map file. This is a tab-delimited file **with a header row**. It must contain columns `mohd_accession`, `r1`, and `r2`. Lines beginning with `#` are treated as comments. Each `mohd_accession` must be unique.

  Example:

  ```
  mohd_accession	r1	r2
  MOHD000001	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
  MOHD000002	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz
  ```

- `input.tmp_dir`: Temporary directory for intermediate files. In some compute environments, `/tmp` is mounted with `noexec`, in which case you should set this to `.tmp`.

#### Resources

```yaml
resources:
  tss_reference: "resources/ENCFF493CCB.bed"
  reference_fasta: "resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
  reference_chromsizes: "resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.chromsizes"
```

These are paths to the reference genome and annotation files required by the pipeline. The TSS reference is used for QC metrics.

#### Parameters

```yaml
params:
  k: 4
  q: 0.05
  p: 0.01
  samtools_exclude_flag: 1804
  samtools_include_flag: 2
```

- `params.k`: Number of alignments to report per read for Bowtie2 (`-k` flag). Default is 4.
- `params.q`: MACS3 q-value threshold for peak calling on the full sample.
- `params.p`: MACS3 p-value threshold for peak calling on pseudoreplicates (used in IDR analysis).
- `params.samtools_exclude_flag`: SAM flag bits to exclude during filtering. The default value of `1804` excludes unmapped reads (4), mate-unmapped reads (8), secondary alignments (256), reads failing QC (512), and PCR duplicates (1024).
- `params.samtools_include_flag`: SAM flag bits to require. The default value of `2` requires proper pairs.

### Pipeline Steps

The pipeline proceeds through the following stages:

#### Preprocessing

- **fastp** — Adapter trimming and quality filtering of raw paired-end FASTQ files. Produces trimmed FASTQs and HTML/JSON QC reports.

#### Alignment

- **bowtie2_build** — Builds a Bowtie2 index from the reference genome. Also generates chromosome sizes and a BED-format chromsizes file used by downstream signal track steps.
- **bowtie2_align** — Aligns trimmed reads to the reference genome with Bowtie2 (`-X 2000 --mm -k 4`). Output is a coordinate-sorted, indexed BAM.

#### Alignment QC

- **pbc_qc** — Computes PCR Bottleneck Coefficients (PBC1, PBC2) and the Non-Redundant Fraction (NRF) from the aligned BAM. These library complexity metrics follow ENCODE standards. Output is a JSON report.

#### Filtering

- **filter** — Marks duplicates with sambamba, applies SAM flag filtering, enforces MAPQ ≥ 30, removes mitochondrial reads, and produces a final filtered BAM. Includes `samtools fixmate` and `samtools quickcheck` for integrity validation.

#### Format Conversion

- **bam_to_tagalign** — Converts the filtered BAM to BEDPE and gzipped tagAlign format for use with MACS3 and pseudoreplicate generation.

#### Signal Track Generation

- **macs3_signal** — Runs MACS3 `callpeak` with `-B` to produce treatment pileup and control lambda bedGraph files.
- **macs3_fold_change_signal** — Computes fold-change signal (FE) over control using `macs3 bdgcmp`, then converts to bigWig.
- **macs3_pvalue_signal** — Computes p-value signal (Poisson) using `macs3 bdgcmp`, then converts to bigWig.
- **bedgraphtobigwig** — Converts the treatment pileup bedGraph to bigWig format using UCSC `bedGraphToBigWig`.

#### Peak Calling

- **macs3_callpeak** — Calls peaks on the full sample tagAlign with MACS3 using the configured q-value threshold. Produces a narrowPeak file.

#### Library QC

- **frag_len** — Computes and plots the fragment length distribution from the filtered BAM. Reports the percentage of fragments in nucleosome-free (< 147 bp), mono-nucleosomal (147–294 bp), di-nucleosomal (294–441 bp), and tri-nucleosomal (> 441 bp) ranges.
- **tss_enrichment** — Calculates TSS enrichment score by aggregating cut-site signal in a ±2 kb window around annotated TSSs. Produces a plot and a JSON metric file.
- **frip_all** — Calculates the Fraction of Reads in Peaks (FRiP) using bedtools intersection of tagAlign reads against called peaks.

#### Reproducibility Analysis

- **pseudoreps** — Splits the BEDPE file into two pseudoreplicates by random shuffling (seed = 61), then converts each to tagAlign format.
- **macs3_pseudoreps** — Calls peaks on each pseudoreplicate tagAlign using MACS3 with the configured p-value threshold.
- **idr** — Runs Irreproducible Discovery Rate (IDR) analysis on pseudoreplicate peak sets. Produces IDR-filtered peaks and a diagnostic plot.
- **overlap_peaks** — Identifies peaks present in both pseudoreplicates using `bedtools intersect`.
- **overlap_peaks_abc** — Identifies peaks from the full sample that overlap each pseudoreplicate peak set by ≥ 50% (reciprocal overlap filtering), following the ABC method.

### Defining Resources

All resource definitions for SLURM rules are located in [profile/slurm/config.yaml](profile/slurm/config.yaml).

```yaml
  bowtie2_align:
    runtime: 720m
    constraint: cascadelake
    slurm_extra: "--exclude=z[1071,1024] --cpu-freq=High-High:Performance"
    slurm_partition: 12hours
    threads: 22
    cpus_per_task: 22
    mem: 64000MB
```

In this instance, we are requesting a job on a SLURM cluster with a runtime of 720 minutes, a specific constraint (i.e. cascadelake), a specific partition (i.e. 12hours), 22 threads, 22 CPUs per task, and 64000 MB of memory.

These resource definitions are specific to the Weng Lab's SLURM cluster. If you are running the pipeline on your own cluster, you will need to adjust these values accordingly. We recommend consulting the [snakemake-executor-plugin-slurm documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for further information.

### Running the Pipeline

With the configuration in place, you can run the pipeline:

```bash
snakemake --workflow-profile profile/slurm
```

> [!NOTE]
> It is recommended that you run the pipeline in dry run mode first (add the `-n` flag). Also note that `snakemake` must be run in the root of this repository.

#### Slack Notifications

The pipeline supports optional Slack notifications. Set the `SLACK_WEBHOOK_URL` environment variable (e.g., via a `.env` file in the project root or by exporting it in your shell):

```bash
export SLACK_WEBHOOK_URL="https://hooks.slack.com/services/..."
```

When set, the pipeline will send start, success, and error notifications to the configured Slack channel.

#### Running Locally

If you cannot run on a SLURM cluster, you can run the pipeline locally. However, this is not a supported use case.

To run locally, the `localrules` directive in [workflow/Snakefile](workflow/Snakefile) already lists all rules. When uncommented (as in the default configuration), these rules will execute on the local machine. Comment out the rules you want to submit to SLURM instead.

### Outputs

The pipeline produces results organized under the `results/` directory:

| Directory | Key Outputs |
|---|---|
| `results/fastp/` | Trimmed FASTQs, HTML and JSON QC reports |
| `results/bowtie2_align/` | Sorted BAMs, flagstat JSON |
| `results/pbc_qc/` | Library complexity JSON (PBC1, PBC2, NRF) |
| `results/filter/` | Filtered BAMs, flagstat JSON |
| `results/bedpe_tagalign/` | BEDPE and tagAlign files |
| `results/macs3_signal/` | Treatment pileup bedGraph, bigWig signal tracks |
| `results/macs3_fold_change_signal/` | Fold-change bigWig signal |
| `results/macs3_pvalue_signal/` | P-value bigWig signal |
| `results/macs3_callpeak/` | narrowPeak files |
| `results/frag_len/` | Fragment length distribution PNG and TXT |
| `results/tss_enrichment/` | TSS enrichment plot PNG and score TXT |
| `results/frip_all/` | FRiP metric TXT |
| `results/pseudoreps/` | Pseudoreplicate tagAlign files |
| `results/macs3_pseudoreps/` | Pseudoreplicate narrowPeak files |
| `results/idr/` | IDR-filtered peaks BED and diagnostic PNG |
| `results/overlap_peaks/` | Overlap peaks BED (intersect method) |
| `results/overlap_peaks_abc/` | Overlap peaks BED (≥ 50% reciprocal overlap) |

### Custom Scripts

The pipeline includes several custom Python scripts located in `workflow/rules/scripts/`:

- **`pbc_qc.py`** — Computes library complexity metrics (PBC1, PBC2, NRF) from a BAM file using `oxbow` and `polars`.
- **`bam_to_tagalign.py`** — Converts a filtered BAM to BEDPE and tagAlign format.
- **`plot_fragment_length_distr.py`** — Extracts fragment lengths from a BAM and generates a distribution plot with nucleosome-region annotations.
- **`tss_enrichment.py`** — Calculates TSS enrichment score by extracting Tn5 cut-site signal around annotated TSSs.
- **`calculate_frip.py`** — Computes the Fraction of Reads in Peaks using bedtools.
- **`bedpe_to_pseudoreps.py`** — Splits a BEDPE file into two pseudoreplicate tagAlign files by random shuffling.
- **`overlap_peaks_abc.py`** — Identifies peaks with ≥ 50% reciprocal overlap across pseudoreplicates using `bioframe`.

## References

> Buenrostro, J., Giresi, P., Zaba, L. et al. _Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position_. Nat Methods 10, 1213–1218, **2013**. https://doi.org/10.1038/nmeth.2688

> Langmead, B., Salzberg, S. _Fast gapped-read alignment with Bowtie 2_. Nat Methods 9, 357–359, **2012**. https://doi.org/10.1038/nmeth.1923

> Zhang, Y., Liu, T., Meyer, C.A. et al. _Model-based Analysis of ChIP-Seq (MACS)_. Genome Biol 9, R137, **2008**. https://doi.org/10.1186/gb-2008-9-9-r137

> Qunhua Li, James B. Brown, Haiyan Huang, Peter J. Bickel. _Measuring reproducibility of high-throughput experiments_. The Annals of Applied Statistics, 5(3), 1752-1779, **2011**.

> Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. _Twelve years of SAMtools and BCFtools_. GigaScience, 10(2), giab008, **2021**. https://doi.org/10.1093/gigascience/giab008

> Artem Tarasov, Albert J. Vilella, Edwin Cuppen, Isaac J. Nijman, Pjotr Prins. _Sambamba: fast processing of NGS alignment formats_. Bioinformatics, 31(12), 2032–2034, **2015**. https://doi.org/10.1093/bioinformatics/btv098

> Chen, Shifu. _fastp 1.0: An Ultra-Fast All-Round Tool for FASTQ Data Quality Control and Preprocessing_. iMeta 4, e70078, **2025**. https://doi.org/10.1002/imt2.70078

> Mölder F, Jablonski KP, Letcher B et al. _Sustainable data analysis with Snakemake_ [version 3; peer review: 2 approved]. F1000Research, 10:33, **2025**. https://doi.org/10.12688/f1000research.29032.3

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
