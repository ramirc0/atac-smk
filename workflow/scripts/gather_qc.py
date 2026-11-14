"""Collect every sample's QC into one table, a row per sample, for the run's report."""

import argparse
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl
from qc_report import (
    count_records,
    read_cutadapt,
    read_fastp,
    read_fingerprint,
    read_flagstat,
    read_picard_metrics,
)


def build_parser():
    """Build the gather_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the gather_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Collect each sample's QC metrics.")
    p.add_argument("--run-dir", required=True, help="Run directory of sample dirs.")
    p.add_argument("--samples", required=True, nargs="+", help="Sample IDs.")
    p.add_argument("--threads", type=int, default=1, help="Samples read at once.")
    p.add_argument("--output", required=True, help="TSV to write.")
    return p


def read_fields(path, column=1):
    """Read a ``name<TAB>value...`` QC file into floats.

    Parameters
    ----------
    path : pathlib.Path
        QC file such as ``frac_mito.qc`` or ``annot_enrich.qc``.
    column : int, default 1
        Field after the name to keep, counting the name as 0.

    Returns
    -------
    dict of str to float
        Value per name.
    """
    fields = (line.split("\t") for line in path.read_text().splitlines())
    return {f[0]: float(f[column]) for f in fields}


def read_nucleosomal(path):
    """Read ENCODE's fragment length QC.

    Parameters
    ----------
    path : pathlib.Path
        ``nucleosomal.qc``: ratio rows hold name, pass, value and message.
        Peak rows hold name, pass and message.

    Returns
    -------
    dict
        ``frac_nfr``, ``nfr_mono_ratio``, and whether the NFR, mono- and
        di-nucleosome peaks are present.
    """
    rows = {r[0]: r for r in (x.split("\t") for x in path.read_text().splitlines())}
    return {
        "frac_nfr": float(rows["Fraction of reads in NFR"][2]),
        "nfr_mono_ratio": float(rows["NFR / mono-nuc reads"][2]),
        "nfr_peak": rows["Presence of NFR peak"][1] == "True",
        "mono_nuc_peak": rows["Presence of Mono-Nuc peak"][1] == "True",
        "di_nuc_peak": rows["Presence of Di-Nuc peak"][1] == "True",
    }


def sample_metrics(run, sample):
    """Collect one sample's metrics, keyed by the ids in ``workflow/report/qc.yaml``.

    Parameters
    ----------
    run : pathlib.Path
        Run directory.
    sample : str
        Sample ID.

    Returns
    -------
    dict
        One row of the metrics table.
    """
    d = run / sample
    fastp = read_fastp(d / f"fastp/{sample}.json")
    aligned = read_flagstat(d / f"bowtie2_align/{sample}.flagstat.json")
    filtered = read_flagstat(d / f"filter/{sample}.flagstat.json")
    pbc = json.loads((d / f"pbc_qc/{sample}.pbc.json").read_text())
    annot = read_fields(d / f"annot_enrich/{sample}.annot_enrich.qc", column=2)
    frip = pl.read_csv(d / f"frip/{sample}.tsv", separator="\t")
    frip = dict(zip(frip["peak_set"], frip["frip"]))
    gc = read_picard_metrics(d / f"gc_bias/{sample}.gc_summary.txt").row(0, named=True)
    return {
        "sample_id": sample,
        "read_pairs": fastp["reads"] // 2,
        "read_length": fastp["read_length"],
        "q30_pct": fastp["q30_pct"],
        "gc_pct": fastp["gc_pct"],
        "duplication_pct": fastp["duplication_pct"],
        "adapter_pct": read_cutadapt(d / f"trim/{sample}.cutadapt.json")["adapter_pct"],
        "alignment_pct": aligned["primary_mapped_pct"],
        "properly_paired_pct": aligned["properly_paired_pct"],
        "frac_mito": read_fields(d / f"frac_mito/{sample}.frac_mito.qc")[
            "frac_mito_reads"
        ],
        "fragments": filtered["primary"] // 2,
        "nrf": pbc["NRF"],
        "pbc1": pbc["PBC1"],
        # ENCODE writes -1 when no location holds exactly two pairs. The portal
        # reads that as infinite.
        "pbc2": math.inf if pbc["PBC2"] == -1 else pbc["PBC2"],
        "tss_enrichment": float(
            (d / f"tss_enrichment/{sample}.tss_enrich.qc").read_text()
        ),
        "frip": frip["overlap"],
        "frip_peaks": frip["peaks"],
        "frip_idr": frip["idr"],
        "frac_dhs": annot["fraction_of_reads_in_universal_DHS_regions"],
        "frac_promoter": annot["fraction_of_reads_in_promoter_regions"],
        "frac_enhancer": annot["fraction_of_reads_in_enhancer_regions"],
        "frac_blacklist": annot["fraction_of_reads_in_blacklist_regions"],
        "synthetic_jsd": read_fingerprint(d / f"jsd/{sample}.jsd.qc")["synthetic_jsd"],
        "at_dropout": gc["AT_DROPOUT"],
        "gc_dropout": gc["GC_DROPOUT"],
        "peaks": count_records(d / f"peaks/{sample}.bfilt.narrowPeak.gz"),
        "overlap_peaks": count_records(
            d / f"overlap/{sample}.overlap.bfilt.narrowPeak.gz"
        ),
        "idr_peaks": count_records(d / f"idr/{sample}.idr.bfilt.narrowPeak.gz"),
        **read_nucleosomal(d / f"frag_len/{sample}.nucleosomal.qc"),
    }


def main(argv=None):
    """Write the run's metrics table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    run = Path(args.run_dir)
    with ThreadPoolExecutor(args.threads) as pool:
        rows = list(pool.map(lambda s: sample_metrics(run, s), args.samples))
    pl.DataFrame(rows).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    main()
