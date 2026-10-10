# CLAUDE.md

ATAC-seq Snakemake workflow on the snakemake-template layout.
README.md holds the layout, steps, outputs, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
pixi run snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml  # ENCODE end-to-end test
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, python-dotenv, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs `default` and `idr`. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV` and `IDR_ENV` resolve `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Outputs must not drift

- Rule env pins match `clarity001/atac-smk:opt`, the image behind the production results. `:latest` is broken: no `bowtie2-build`, and `import polars` fails.
- Tag `pre-template` is the last commit before the template port. Its outputs and this workflow's are byte-identical.
- README.md "Known issues" lists legacy bugs kept on purpose. A fix changes outputs and needs the user's approval.

## Gotchas

- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Plots use poikilos `plain`. `pk.save_figure()` writes SVG and PNG. Rules MUST declare both outputs.
- The `default` env sets `MPLBACKEND=Agg` and `LC_ALL=C`. The `idr` env sets neither, like the conda env it replaced.
- conda's `bowtie2-build` lacks SAIS. `--large-index` keeps the GRCh38 index as `.bt2l`.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- CLI `--set-threads` replaces the profile's whole mapping. Pass every rule.

## Verifying a change

No test suite. Diff against the legacy pipeline on a real sample:

1. Take the first 3M read pairs of a sample's raw FASTQs.
2. Run `git archive pre-template` with the `:opt` image and its conda IDR env.
3. Run this workflow with the same per-rule threads.
4. Run `fastp` and `bowtie2_align` single-threaded in both. Multithreaded bowtie2 reorders reads, and sambamba then keeps different duplicates.
5. Compare BAM records with `samtools view`, gzipped files decompressed, and MACS3 `.xls` files without `#` lines.
