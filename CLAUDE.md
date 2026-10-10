# CLAUDE.md

ATAC-seq Snakemake workflow on the snakemake-template layout.
README.md holds the layout, steps, outputs, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local
pixi run test                                # chr21 fixture vs ENCODE's task scripts (TMPDIR off /tmp)
pixi run snakemake --profile profiles/slurm --configfile .tests/integration/config.yaml  # ENCODE end-to-end test
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, python-dotenv, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs `default`, `idr` and `test`. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`. `pixi lock` keeps locked versions; `pixi update` takes the newest.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV` and `IDR_ENV` resolve `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Match ENCODE

- Every step reproduces ENCODE ATAC-seq pipeline v2.2.3 (`~/Git_Repos/atac-seq-pipeline`), single-replicate flow. README.md "Differences from ENCODE" lists the deliberate gaps.
- `.tests/fixture/test_encode.py` checks each step against ENCODE's own task script on the same inputs. A step change MUST keep its test passing.
- A bug fix starts with a test that fails on the old code.
- Correctness beats parity. Where ENCODE's tool has a known bug, use the fixed version and list the gap in README.md.
- Tag `pre-template` is the legacy pipeline. Tag `post-template` is the same outputs on the template layout, the last commit before the ENCODE port.

## Gotchas

- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- In `-k` mode bowtie2 gives uniquely aligned pairs MAPQ 255, and multimappers ordinary MAPQs.
- sambamba markdup's hash and overflow sizes change which duplicate survives. Keep ENCODE's.
- pytest's diff of long line lists takes minutes. Compare them with `first_difference`.
- Plots use poikilos `plain`. `pk.save_figure()` writes SVG and PNG. Rules MUST declare both outputs.
- The `default` env sets `MPLBACKEND=Agg` and `LC_ALL=C`. The `idr` env sets `LC_ALL=C`, because its rule sorts.
- conda's `bowtie2-build` lacks SAIS. `--large-index` keeps the GRCh38 index as `.bt2l`.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- CLI `--set-threads` replaces the profile's whole mapping. Pass every rule.

## Verifying a change

1. `pixi run lint` and `pixi run test`. The login node's `/tmp` fills up: set `TMPDIR` to `/data` scratch.
2. For a big change, run the integration test on SLURM, then `.tests/integration/compare_encode.py`.
