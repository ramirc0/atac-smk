# atac-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for paired-end
ATAC-seq. It reproduces the single-replicate flow of the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
v2.2.3, step by step. Reads are trimmed and aligned. The filtered alignments
feed signal tracks, peak calls, IDR and naive overlap on pseudoreplicates, and
ENCODE's library and enrichment QC. Fold change is z-scored over the ENCODE
cCRE anchors, as the ENCODE4 cCRE pipeline does.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, r1, r2
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, path constants
    preprocess.smk       trimming, bowtie2 indexes and alignment, filtering, tagAlign
    peaks.smk            pseudoreplicates, MACS3 peaks, blacklist filter, IDR, overlap
    signal.smk           fold-change, p-value and pileup bigWigs, anchor z-scores
    qc.smk               PBC/NRF, fragment lengths, TSS, FRiP, mito, annotation, GC, JSD
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule envs (pixi workspace): default, idr and test
.tests/
  integration/           chr21 fixture run, step by step against ENCODE's task scripts
  validation/            four ENCODE replicates on GRCh38, compared with ENCODE's QC
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env, lint and test tasks
pyproject.toml           ruff + numpydoc config
```

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
cp config/samples.example.tsv config/samples.tsv   # then edit
pixi shell                                         # launcher env
snakemake -n -p --profile profiles/local           # dry run
snakemake --profile profiles/local                 # local
snakemake --profile profiles/slurm                 # SLURM
```

Snakemake MUST run from the repository root. The SLURM profile holds Weng Lab
cluster settings, such as partitions and the `cascadelake` constraint. Adapt them
for another cluster with the
[SLURM executor docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).
Tune `set-resources` from the `benchmarks/<run_id>/` TSVs after a first run.

## Testing

Both tests use data on the Weng Lab cluster, so they run only there.

`pixi run test` runs the integration tests. The fixture holds the chr21 and chrM
reads of ENCSR095QNB replicates 1 and 2, named `S1` and `S10`. The tests run the
workflow on it. Then they run ENCODE's own task scripts, from the
`encodedcc/atac-seq-pipeline:v2.2.3` image, on the same inputs, and compare.
The z-score test runs the ENCODE4 cCRE pipeline's script in that image, from a
clone of [ENCODE-cCREs](https://github.com/weng-lab/ENCODE-cCREs). Set these up
once:

```bash
pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
apptainer pull /data/zusers/ramirezc/atac-smk-refs/containers/atac-seq-pipeline_v2.2.3.sif \
    docker://encodedcc/atac-seq-pipeline:v2.2.3   # or set ENCODE_SIF
git clone https://github.com/weng-lab/ENCODE-cCREs ~/Git_Repos/ENCODE-cCREs
pixi run test
```

The validation runs four ENCODE replicates on full GRCh38, then prints
their QC beside ENCODE's published values:

```bash
snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml
pixi run --manifest-path workflow/envs/pixi.toml python .tests/validation/compare_encode.py
```

## Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its raw paired-end FASTQs. Lines that start with `#` are comments.
Each `sample_id` MUST be unique.

The references are ENCODE's GRCh38 files. Download them from
`https://www.encodeproject.org/files/<accession>/@@download/<accession>.bed.gz`.

| Key | File |
|---|---|
| `references.fasta` | GRCh38 no-alt analysis set FASTA (`GRCh38_no_alt_analysis_set_GCA_000001405.15`) |
| `references.tss` | TSS, ENCFF493CCB |
| `references.blacklist` | Exclusion list v2, ENCFF356LFX |
| `references.dhs` | Universal DHS regions, ENCFF304XEX |
| `references.promoters` | Promoters, ENCFF140XLU |
| `references.enhancers` | Enhancers, ENCFF212UAV |
| `references.anchors` | cCRE registry V4 anchors (rDHSs and TF rClusters), `GRCh38-Anchors.bed` |

The anchors are not on the portal. The Weng Lab keeps them at
`/data/projects/encode/Registry/V4/GRCh38/GRCh38-Anchors.bed`: 2,950,228 rows,
MD5 `9c5d57a3de8f09b27e9c00661801eddb`.

## Configuration

| Key | Meaning |
|---|---|
| `samples` | Sample sheet path. |
| `outdir` | Results root. |
| `run_id` | Run name. Outputs go under `results/<run_id>/`. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` that most rules use. |
| `mito_chr_name` | Mitochondrial contig. |
| `tmp_dir` | Scratch for sorts and per-job temp directories. Use `.tmp` where `/tmp` is small. |
| `references.*` | See Inputs. |
| `params.multimapping` | Alignments a read may have. bowtie2 runs with `-k` one higher. |
| `params.pval_thresh` | MACS3 p-value cutoff. |
| `params.cap_num_peak` | Peaks kept per call, by p-value. |
| `params.smooth_win` | MACS3 `--extsize`. `--shift` is minus half of it. |
| `params.idr_thresh` | IDR threshold on global IDR. |

## Steps

Preprocessing, in `preprocess.smk`:

- `fastp_qc` writes fastp's HTML and JSON report on the raw reads. It trims nothing.
- `trim` detects each mate's adapter as ENCODE does and trims it with
  `cutadapt -e 0.1 -m 5`.
- `bowtie2_build` and `bowtie2_build_mito` index the genome and its mitochondrial contig.
- `bowtie2_align` aligns with `-k 5 -X 2000 --mm` into a sorted, indexed BAM.
- `filter` keeps proper pairs. Read names with more than 4 alignments are dropped,
  others keep their primary one at any MAPQ. Then sambamba marks duplicates,
  which are removed with the mitochondrial reads.
- `tagalign` writes one Tn5-shifted line per mate, in name order.

Peaks, in `peaks.smk`:

- `pseudoreps` shuffles whole read pairs with ENCODE's seeded keystream and halves them.
- `macs3` calls the sample's peaks and writes the bedGraphs for signal.
  `macs3_pseudorep` calls each pseudoreplicate's peaks. Both use
  `-p 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits`,
  sort by p-value and keep 300,000 peaks.
- `blacklist_filter` drops peaks that touch the exclusion list or lie off chr1-22, X and Y.
- `idr` scores the pseudoreplicate peaks against the sample's own. It keeps
  peaks at global IDR ≤ 0.05. `plot_idr` draws the IDR panels.
- `overlap` keeps sample peaks that overlap a peak in each pseudoreplicate by
  half of either peak.

Signal, in `signal.smk`:

- `signal` writes fold-change and p-value bigWigs with `macs3 bdgcmp`, plus a
  pileup bigWig in signal per million reads.
- `zscore` averages fold change over each anchor with bigtools
  `bigwigaverageoverbed`. It z-scores log10 of the non-zero averages with their
  mean and population standard deviation, as the ENCODE4 cCRE pipeline does.
  Anchors without signal get -10.

QC, in `qc.smk`:

- `pbc_qc` computes NRF, PBC1 and PBC2 on the filtered reads before duplicate removal.
- `frag_len` runs Picard CollectInsertSizeMetrics and ENCODE's NFR and
  nucleosome checks.
- `tss_enrichment` ports ENCODE's metaseq TSS enrichment.
- `frip` counts tagAlign reads in each filtered peak set.
- `frac_mito` aligns every read to the mitochondrial contig alone, as ENCODE does.
- `annot_enrich` counts reads in DHS, promoter, enhancer and blacklist regions.
- `gc_bias` runs Picard CollectGcBiasMetrics. `jsd` runs deeptools plotFingerprint.

## Outputs

Each sample writes to `results/<run_id>/<sample>/<step>/`. The indexes and
chromosome sizes do not depend on the run. They go to `results/refs/`. Logs and
benchmarks mirror the run layout under `logs/<run_id>/` and `benchmarks/<run_id>/`.
Trimmed reads, bedGraphs and unfiltered IDR and overlap peaks are deleted once used.

| Step directory | Key outputs |
|---|---|
| `fastp/` | fastp HTML and JSON report |
| `trim/` | cutadapt JSON report |
| `bowtie2_align/` | Sorted BAM, flagstat JSON |
| `filter/` | Filtered BAM and index, flagstat JSON |
| `tagalign/` | Tn5-shifted tagAlign |
| `pseudoreps/` | Pseudoreplicate tagAligns |
| `peaks/` | Sample and pseudoreplicate narrowPeaks, filtered sample peaks |
| `idr/` | IDR peaks (filtered), all IDR scores, IDR plot |
| `overlap/` | Naive overlap peaks (filtered) |
| `signal/` | Fold-change, p-value and pileup bigWigs |
| `zscore/` | Per-anchor fold-change z-score, mean signal and rank, as parquet |
| `pbc_qc/` | NRF, PBC1 and PBC2 JSON |
| `frag_len/` | ENCODE's `nucleosomal.qc`, distribution plot |
| `tss_enrichment/` | ENCODE's `tss_enrich.qc`, profile plot |
| `frip/` | FRiP per peak set |
| `frac_mito/` | Mitochondrial read fraction |
| `annot_enrich/` | Read fractions per annotation |
| `gc_bias/` | Picard GC metrics, plot |
| `jsd/` | deeptools fingerprint metrics, plot |

## Slack notifications

Set `SLACK_WEBHOOK_URL` in the shell or in a `.env` file in the repository root.
The pipeline then posts a Slack message when a run starts and another when it
ends.

```bash
export SLACK_WEBHOOK_URL="https://hooks.slack.com/services/..."
```

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main and its plugins
for pixi deployment and SLURM. The rule envs live in a second pixi workspace,
`workflow/envs/`. Its `pixi.lock` pins every package.
`default` runs most rules. `idr` runs IDR on the Python 3.9 / NumPy 1.19 stack it
needs. `test` adds pytest and seqkit. pixi MUST be on `PATH`. An `onstart` hook
installs the envs before the first job. Parallel SLURM jobs would otherwise race
to install them.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). `pixi update` moves the lock to
the newest versions the manifest allows. Any change to the manifest or lock
reruns every rule. So does moving the workdir, because the env hash includes the
workspace's absolute path.

## Conventions

The workflow follows the
[Nextstrain Snakemake style guide](https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html)
and [Snakemake's best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html):

- No `run:` blocks. Custom logic lives in `workflow/scripts/`, called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation uses `:q` quoting.
- Config reaches `shell` through `params:` (a `lambda`), never directly.
- `input:` and `output:` are literal path strings. All paths are relative.
- Rules tee logs to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has `log:` and `benchmark:`. Each names its pixi env in `software:`.
- Plots use [poikilos](https://github.com/ramirc0/poikilos) with the `plain`
  theme. `pk.save_figure()` writes the SVG and PNG pair. Rules declare both.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and
  `numpydoc lint`. Neither tool parses `.smk` files.

## Differences from ENCODE

- ENCODE's `assign_multimappers.py` drops the last read name of every BAM. This does not.
- MACS3 3.0.5 replaces MACS2 2.2.4. Peak regions match. Summits differ where
  MACS3 fixed summit selection bugs.
- sambamba marks duplicates, with ENCODE's settings. ENCODE's default is Picard.
- Fragment length QC reads every record. ENCODE stops after 5M, mostly chr1.
- One MACS3 call per sample gives peaks and signal. ENCODE runs the same call twice.
- There is no replicate pooling or true-replicate IDR. Each sample is one replicate.
- Anchor z-scores use bigtools' unrounded means. The ENCODE4 cCRE pipeline's UCSC
  tool rounds them to six significant digits. Z-scores agree to 1e-5, and anchors
  tied in its ranks may split.

## References

Buenrostro, Jason D., Paul G. Giresi, Lisa C. Zaba, Howard Y. Chang, and William J. Greenleaf. 2013. “Transposition of Native Chromatin for Fast and Sensitive Epigenomic Profiling of Open Chromatin, DNA-Binding Proteins and Nucleosome Position.” Nature Methods 10 (12): 1213–18. https://doi.org/10.1038/nmeth.2688.

Chen, Shifu. 2025. “Fastp 1.0: An Ultra-Fast All-Round Tool for FASTQ Data Quality Control and Preprocessing.” iMeta 4 (5): e70078. https://doi.org/10.1002/imt2.70078.

Danecek, Petr, James K. Bonfield, Jennifer Liddle, et al. 2021. “Twelve Years of SAMtools and BCFtools.” GigaScience 10 (2): giab008. https://doi.org/10.1093/gigascience/giab008.

ENCODE Project Consortium, Jill E. Moore, Michael J. Purcaro, et al. 2020. “Expanded Encyclopaedias of DNA Elements in the Human and Mouse Genomes.” Nature 583 (7818): 699–710. https://doi.org/10.1038/s41586-020-2493-4.

Huey, Jack D., and Nezar Abdennur. 2024. “Bigtools: A High-Performance BigWig and BigBed Library in Rust.” Bioinformatics (Oxford, England) 40 (6): btae350. https://doi.org/10.1093/bioinformatics/btae350.

Langmead, Ben, and Steven L. Salzberg. 2012. “Fast Gapped-Read Alignment with Bowtie 2.” Nature Methods 9 (4): 357–59. https://doi.org/10.1038/nmeth.1923.

Li, Qunhua, James B. Brown, Haiyan Huang, and Peter J. Bickel. 2011. “Measuring Reproducibility of High-Throughput Experiments.” The Annals of Applied Statistics 5 (3): 1752–79. https://doi.org/10.1214/11-AOAS466.

Martin, Marcel. 2011. “Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads.” EMBnet.Journal 17 (1): 10–12. https://doi.org/10.14806/ej.17.1.200.

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, et al. 2021. “Sustainable Data Analysis with Snakemake.” F1000Research 10 (January): 33. https://doi.org/10.12688/f1000research.29032.3.

Ramírez, Fidel, Devon P. Ryan, Björn Grüning, et al. 2016. “deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis.” Nucleic Acids Research 44 (W1): W160–65. https://doi.org/10.1093/nar/gkw257.

Tarasov, Artem, Albert J. Vilella, Edwin Cuppen, Isaac J. Nijman, and Pjotr Prins. 2015. “Sambamba: Fast Processing of NGS Alignment Formats.” Bioinformatics (Oxford, England) 31 (12): 2032–34. https://doi.org/10.1093/bioinformatics/btv098.

Zhang, Yong, Tao Liu, Clifford A. Meyer, et al. 2008. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9 (9): R137. https://doi.org/10.1186/gb-2008-9-9-r137.

## Questions

If you have any questions or would like to provide constructive feedback, please open an issue or reach out to the MOHD DACC.
