from pathlib import Path
fromt dotenv import load_dotenv
import polars as pl
import math
import glob
import sys
import os

load_dotenv()

Path(config["input"]["tmp_dir"]).mkdir(parents=True, exist_ok=True)
sample_map_df = pl.read_csv(config["input"]["sample_map"], separator="\t", comment_prefix="#")
SAMPLES = sample_map_df["mohd_accession"].to_list()
read1 = sample_map_df["r1"].to_list()
read2 = sample_map_df["r2"].to_list()
FIRST_SAMPLE = read1[0]
READ_LOOKUP = dict(zip(SAMPLES, zip(read1, read2)))

def get_fastqs(wildcards):
    r1, r2 = READ_LOOKUP[wildcards.sample]
    return {"r1": r1, "r2": r2}

ATTEMPT_MODIFIERS = {
    1: {"slurm_partition": "4hours", "runtime": 240, "mem_mb_multiplier": 1.0},
    2: {"slurm_partition": "12hours", "runtime": 720, "mem_mb_multiplier": 1.5},
    3: {"slurm_partition": "5days", "runtime": 7200, "mem_mb_multiplier": 2.0},
}

def get_rsources(wildcards, attempt, input):
    attempt_modifiers = ATTEMPT_MODIFIERS[attempt]
    mem = int(attempt_modifiers["mem_mb_multiplier"]) * input.resources.mem_mb
    return {
        "slurm_partition": attempt_modifiers["slurm_partition"],
        "runtime": attempt_modifiers["runtime"],
        "mem": f"{mem}MB",
    }

def get_macs3_pseudoreps_signal(wildcards):
    return expand(rules.macs3_pseudoreps.output.signal, sample=[wildcards.sample],pseudorep=["pseudorep1","pseudorep2"])
