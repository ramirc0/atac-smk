#!/usr/bin/env python3
import re
import sys

def extract_failed_job_ids(filepath):
    """Extract all SLURM job IDs that failed from a Snakemake log file."""
    failed_jobs = []
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    pattern = r"SLURM-job '(\d+)' failed"
    
    matches = re.findall(pattern, content)
    failed_jobs.extend(matches)
    
    return failed_jobs

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <logfile>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    
    try:
        job_ids = extract_failed_job_ids(filepath)
        
        if job_ids:
            print("Failed SLURM job IDs:")
            print(",".join(job_ids))
        else:
            print("No failed jobs found")
            
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found")
        sys.exit(1)