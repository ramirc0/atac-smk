import subprocess
import argparse
import glob
import os
import yaml
from collections import defaultdict
from pathlib import Path
import json

def _parse_flagstat(flagstat_file: str) -> int:
    """Extract mapped reads from flagstat file"""
    with open(flagstat_file) as f:
        for line in f:
            if "mapped (" in line:
                return int(line.split()[0])

def generate_flagstat_yaml(outdir: str, custom_data_filename: str):
    """Generate MultiQC custom content YAML"""
    
    flagstat_dirs = [os.path.join(outdir, "bowtie2_align"), os.path.join(outdir, "picard")]
    multiqc_dir = os.path.join(outdir, "multiqc")
    
    data = defaultdict(dict)
    flagstat_files = []
    
    for dir in flagstat_dirs:
        flagstat_files.extend(glob.glob(os.path.join(dir, "*.flagstat")))
    
    for flagstat_file in flagstat_files:
        sample, _, _ = os.path.basename(flagstat_file).split(".")
        sample = sample.split("-")[0]
        step = Path(flagstat_file).parent.name
        
        filter_mapped = _parse_flagstat(flagstat_file)
        if filter_mapped:
            data[sample][f"reads_mapped_after_{step}"] = filter_mapped / 1_000_000
        
    output = {
        "id": "user_defined_1",
        "description": "Addition Flagstat Reports included by the user",
        "plot_type": "generalstats",
        "headers": {
            "reads_mapped_after_filter": {
                "title": "Reads mapped (samtools filter)",
                "description": "Number of mapped reads after samtools filtering step",
                "format": "{:.1f} M",
            },
            "reads_mapped_after_picard": {
                "title": "Reads mapped (picard)",
                "description": "Number of mapped reads after Picard deduplication",
                "format": "{:.1f} M",
            },
        },
        "data": dict(sorted(data.items()))
    }
    
    with open(os.path.join(multiqc_dir, custom_data_filename), "w") as f:
        yaml.dump(output, f, default_flow_style=False, sort_keys=False)
        
def generate_fastp_yaml(outdir: str, custom_data_filename: str):
    """Generate MultiQC custom content YAML"""
    
    fastp_dir = f"{outdir}/fastp"
    multiqc_dir = f"{outdir}/multiqc"
    
    data = defaultdict(dict)
    
    for fastp_file in glob.glob(os.path.join(fastp_dir, "*.json")):
        sample = Path(fastp_file).stem
        fastp_json = json.load(open(fastp_file))
        data[sample]["reads_before_filtering"] = fastp_json["summary"]["before_filtering"]["total_reads"] / 1_000_000
    
    output = {
        "id": "user_defined_2",
        "description": "Additional fastp fields included by the user",
        "plot_type": "generalstats",
        "headers": {
            "reads_before_filtering": {
                "title": "Reads Before Filtering",
                "description": "Total reads before filtering (millions)",
                "format": "{:.1f} M",
            },
        },
        "data": dict(sorted(data.items()))
    }
    
    with open(os.path.join(multiqc_dir, custom_data_filename), "w") as f:
        yaml.dump(output, f, default_flow_style=False, sort_keys=False)
    
def generate_file_list(macs3_qvalue: float, outdir:str, file_list_filename: str):
    """Generate file list for MultiQC"""
    multiqc_dir = os.path.join(outdir, "multiqc")
    file_list = os.path.join(multiqc_dir, file_list_filename)
    
    with open(file_list, "w") as f:
        analysis_dir = {
            f'results/macs3_callpeak/{macs3_qvalue}': '*.xls',
            'results/bowtie2_align': "*.flagstat",
            'results/picard': "*.txt",
            'results/fastp': "*.json",
            'results/multiqc': "*_mqc.yaml"
        }
        for dir, pattern in analysis_dir.items():
            files = glob.glob(os.path.join(dir, pattern), recursive=True)
            for file in files:
                f.write(file + "\n")
    
def run_multiqc(multiqc_config: str, outdir: str, file_list_filename: str, overwrite: bool):
    """Run MultiQC"""
    multiqc_dir = os.path.join(outdir, "multiqc")
    file_list = os.path.join(multiqc_dir, file_list_filename)
    cmd = ["multiqc", "--config", multiqc_config, "--file-list", file_list, "--outdir", multiqc_dir]
    if overwrite:
        cmd.append("--force")
    subprocess.check_call(cmd)

def main(args):
    outdir = args.outdir
    multiqc_config = args.multiqc_config
    macs3_qvalue = args.macs3_qvalue
    flagstat_mqc_yaml = args.flagstat_mqc_yaml
    fastp_mqc_yaml = args.fastp_mqc_yaml
    file_list = args.file_list
    overwrite = args.overwrite
    
    multiqc_dir = os.path.join(outdir, "multiqc")
    os.makedirs(multiqc_dir, exist_ok=True)
    
    generate_flagstat_yaml(outdir=outdir, custom_data_filename=flagstat_mqc_yaml)
    generate_fastp_yaml(outdir=outdir, custom_data_filename=fastp_mqc_yaml)
    generate_file_list(macs3_qvalue=macs3_qvalue, outdir=outdir, file_list_filename=file_list)
    run_multiqc(multiqc_config=multiqc_config, outdir=outdir, file_list_filename=file_list, overwrite=overwrite)

if __name__ == "__main__":
    # multiqc {params.multiqc_conf} {RESULTS_DIR} --dirs-depth 1 --outdir {RESULTS_DIR}/multiqc
    parser = argparse.ArgumentParser(description="Prepare custom QC data (*_mqc.yaml) and file list for MultiQC.")
    parser.add_argument("--multiqc_config", help="Path to the MultiQC config file")
    parser.add_argument("--outdir", help="Path to the output directory")
    parser.add_argument("--macs3_qvalue", type=float, help="MACS3 qvalue threshold. Determines which set of peaks will be included in the MultiQC report.")
    parser.add_argument("--flagstat_mqc_yaml", required=False, default="flagstat_mqc.yaml", help="Name of the flagstat YAML file")
    parser.add_argument("--fastp_mqc_yaml", required=False, default="fastp_mqc.yaml", help="Name of the fastp YAML file")
    parser.add_argument("--file_list", required=False, default="file_list.txt", help="Name of the file list file")
    parser.add_argument("--overwrite", default=False, action="store_true", help="Overwrite any existing multiqc data")
    args = parser.parse_args()
    main(args)