import polars as pl
import os
import glob
import sys

metadata = pl.read_csv(config['input_files']["metadata"], separator="\t")

GENOMES = list(config['genomes'].keys())
SAMPLES = metadata["Sample"].to_list()
QVALS = config['qvals']
MIN_QVAL = min(QVALS)
LOGDIR = config['logdir']
read1 = metadata["R1"].to_list()
FIRST_SAMPLE = read1[0]
read2 = metadata["R2"].to_list()

READ_LOOKUP = dict(zip(SAMPLES, zip(read1, read2)))

RESULTS_DIR = config['results_dir']

ATTEMPT_MODIFIERS = {
    1: {"slurm_partition": "4hours", "runtime": 240, "mem_mb_multiplier": 1.0},
    2: {"slurm_partition": "12hours", "runtime": 720, "mem_mb_multiplier": 1.5},
    3: {"slurm_partition": "5days", "runtime": 7200, "mem_mb_multiplier": 2.0},
}

def get_slurm_partition(wildcards, attempt):
    return ATTEMPT_TO_QUEUE[attempt]["slurm_partition"]

def get_runtime(wildcards, attempt):
    return ATTEMPT_TO_QUEUE[attempt]["runtime"]

def get_memory(wildcards, attempt, input):
    return int(ATTEMPT_TO_QUEUE[attempt]["mem_mb_multiplier"] * input.size_mb)