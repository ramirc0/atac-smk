# atac-smk

A [Snakemake](https://github.com/snakemake/snakemake) pipeline for paired-end
ATAC-seq. It follows the [ENCODE ATAC-seq processing standards](https://www.encodeproject.org/atac-seq/).
Reads are trimmed and aligned. The filtered alignments feed signal tracks, peak
calls, library QC and IDR on pseudoreplicates.

## Layout

```
config/
  config.yaml.template   parameters + paths; copy to config.yaml (gitignored)
  samples.example.tsv    one row per sample: sample_id, r1, r2
workflow/
  Snakefile              includes the rule files, defines rule all
  rules/
    common.smk           config, sample sheet, path constants
    preprocess.smk       fastp, bowtie2 index and alignment, filtering, tagAlign
    signal.smk           MACS3 pileup, fold-change and p-value bigWigs
    peaks.smk            MACS3 peaks, pseudoreplicates, IDR, overlap peaks
    qc.smk               PBC/NRF, fragment lengths, TSS enrichment, FRiP
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule envs (pixi workspace): default and idr
.tests/integration/      ENCODE end-to-end test: config + sample sheet
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + lint task
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

## Integration test

`.tests/integration/` runs four ENCODE ATAC-seq replicates end to end on full
GRCh38. Its inputs live on the Weng Lab cluster, so it runs only there.

```bash
snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml
```

Outputs go under `results/encode-integration/`. Compare the QC with ENCODE's
for the same replicates.

## Inputs

The sample sheet is a TSV with a header row. `sample_id` names the sample. `r1`
and `r2` are its raw paired-end FASTQs. Lines that start with `#` are comments.
Each `sample_id` MUST be unique.

The references live in `resources/` (gitignored):

- the GRCh38 FASTA, no-alt analysis set
- a TSS BED, e.g. `ENCFF493CCB.bed`

## Configuration

| Key | Meaning |
|---|---|
| `samples` | Sample sheet path. |
| `outdir` | Results root. |
| `run_id` | Run name. Outputs go under `results/<run_id>/`. |
| `pixi_env` | Environment of `workflow/envs/pixi.toml` that most rules use. |
| `tmp_dir` | Scratch for `samtools sort` and IDR. Use `.tmp` where `/tmp` is mounted `noexec`. |
| `references.fasta` | Reference genome FASTA. |
| `references.tss` | TSS BED for TSS enrichment. |
| `params.k` | Alignments bowtie2 reports per read (`-k`). |
| `params.q` | MACS3 q-value cutoff for full-sample peaks. |
| `params.p` | MACS3 p-value cutoff for pseudoreplicate peaks. |
| `params.samtools_exclude_flag` | SAM flags to drop. `1804` drops unmapped (4), mate unmapped (8), secondary (256), QC fail (512) and duplicate (1024) reads. |
| `params.samtools_include_flag` | SAM flags to require. `2` keeps proper pairs. |

## Steps

Preprocessing, in `preprocess.smk`:

- `fastp` trims adapters and low-quality bases. It writes HTML and JSON reports.
- `bowtie2_build` indexes the reference and writes its chromosome sizes.
- `bowtie2_align` aligns with `-X 2000 --mm -k 4` into a sorted, indexed BAM.
- `filter` marks duplicates with sambamba and keeps proper pairs with MAPQ ≥ 30.
  It drops every SAM line that mentions `chrM`.
- `bam_to_tagalign` writes a BEDPE and a Tn5-shifted, sorted tagAlign.

Signal tracks, in `signal.smk`:

- `macs3_signal` runs `macs3 callpeak -B` for the pileup and control lambda.
- `bedgraphtobigwig` converts the pileup to bigWig.
- `macs3_fold_change_signal` and `macs3_pvalue_signal` run `macs3 bdgcmp`
  (`-m FE` and `-m ppois`) and convert the result to bigWig.

Peaks, in `peaks.smk`:

- `macs3_callpeak` calls full-sample peaks at the `params.q` cutoff.
- `pseudoreps` shuffles the BEDPE rows (seed 61) and splits them into two
  pseudoreplicate tagAligns.
- `macs3_pseudoreps` calls peaks on each pseudoreplicate at the `params.p` cutoff.
- `idr` scores the pseudoreplicate peaks with IDR and plots the result.
- `overlap_peaks` keeps pseudoreplicate 1 peaks that overlap a pseudoreplicate 2
  peak.
- `overlap_peaks_abc` keeps full-sample peaks that overlap a peak in each
  pseudoreplicate. The overlap covers at least half of either peak.

QC, in `qc.smk`:

- `pbc_qc` computes PBC1, PBC2 and NRF from the unfiltered BAM.
- `frag_len` counts fragment lengths and plots them. It reports the share of
  nucleosome-free, mono-, di- and tri-nucleosomal fragments.
- `tss_enrichment` scores Tn5 cut sites in a ±2 kb window around each TSS.
- `frip_all` computes the fraction of tagAlign reads in peaks.

## Outputs

Each sample writes to `results/<run_id>/<sample>/<step>/`. The bowtie2 index and
chromosome sizes do not depend on the run. They go to `results/refs/bowtie2_build/`.
Logs and benchmarks mirror the run layout under `logs/<run_id>/` and
`benchmarks/<run_id>/`.

| Step directory | Key outputs |
|---|---|
| `fastp/` | Trimmed FASTQs, HTML and JSON reports |
| `bowtie2_align/` | Sorted BAM, flagstat JSON |
| `pbc_qc/` | PBC1, PBC2 and NRF JSON |
| `filter/` | Filtered BAM, flagstat JSON |
| `bedpe_tagalign/` | BEDPE, gzipped tagAlign |
| `macs3_signal/` | Pileup and control lambda bedGraphs, pileup bigWig |
| `macs3_fold_change_signal/` | Fold-change bedGraph and bigWig |
| `macs3_pvalue_signal/` | P-value bedGraph and bigWig |
| `macs3_callpeak/` | narrowPeak |
| `frag_len/` | Distribution plot (SVG, PNG), range counts TXT |
| `tss_enrichment/` | Profile plot (SVG, PNG), score JSON in a TXT |
| `frip_all/` | FRiP TXT |
| `pseudoreps/` | Pseudoreplicate tagAligns |
| `macs3_pseudoreps/` | Pseudoreplicate narrowPeaks |
| `idr/` | IDR BED, IDR plot PNG |
| `overlap_peaks/` | Overlap peaks BED |
| `overlap_peaks_abc/` | Overlap peaks BED (50% rule) |

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
needs. pixi MUST be on `PATH`. Jobs install the envs on first use.

The tool versions match `clarity001/atac-smk:opt`, the image that produced the
production results.

Add a dependency with `pixi add <pkg>` in `workflow/envs/`. After a hand edit,
run `pixi lock` there (rules use `locked=True`). Any change to the manifest or
lock reruns every rule. So does moving the workdir, because the env hash includes
the workspace's absolute path.

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

## Known issues

These legacy bugs stay to keep outputs identical to the production results.
A fix changes outputs.

- `bedpe_to_pseudoreps.py` builds pseudoreplicate 1's mate-2 reads from
  pseudoreplicate 2's rows. Pseudoreplicate 1's own mate-2 reads are dropped.
  On chr21 of one sample, half of pseudoreplicate 2 also appears in
  pseudoreplicate 1. The overlap inflates IDR and overlap-peak reproducibility.
- `bam_to_tagalign.py` prints errors and carries on. A failed step can leave
  truncated outputs in a job that reports success.
- `calculate_frip.py` pipes `bedtools intersect` into `wc -l` without
  `pipefail`. A bedtools failure gives a FRiP of 0.
- `idr` runs `rm <tmp_dir>/<sample>*`. It deletes other samples' scratch files
  when sample IDs prefix one another, such as `S1` and `S10`.
- `plot_fragment_length_distr.py` uses inclusive ranges. Fragments of 147, 294
  and 441 bp count in two classes, so the shares can sum past 100%. Fragments
  over 1000 bp count only in the total.
- `tss_enrichment.py` adds a read's cut site once per pileup column the read
  covers. Each read is weighted by its aligned length in the window.
- `pbc_qc.py` takes `pos + len(seq)` as the read end and ignores the CIGAR.

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
