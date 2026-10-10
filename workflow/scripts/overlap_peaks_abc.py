"""Keep full-sample peaks that overlap a peak in each pseudoreplicate by at least 50%."""

import argparse

import bioframe as bf
import pandas as pd

PEAK_COLS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
]


def build_parser():
    """Build the overlap_peaks_abc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the overlap_peaks_abc.py flags.
    """
    p = argparse.ArgumentParser(
        description="Find peaks that overlap both pseudoreplicates by at least 50%."
    )
    p.add_argument("--peaks", required=True, help="Full-sample narrowPeak.")
    p.add_argument(
        "--pseudorep1", required=True, help="First pseudoreplicate narrowPeak."
    )
    p.add_argument(
        "--pseudorep2", required=True, help="Second pseudoreplicate narrowPeak."
    )
    p.add_argument("--output", required=True, help="Overlapping peaks BED to write.")
    return p


def filter_50percent_overlap(df1, df2):
    """Keep the rows of ``df1`` that overlap a row of ``df2`` by at least 50%.

    Either interval may supply the 50%: the overlap must cover half of the
    ``df1`` peak or half of the ``df2`` peak.

    Parameters
    ----------
    df1 : pandas.DataFrame
        Peaks to filter, with the narrowPeak columns.
    df2 : pandas.DataFrame
        Peaks to overlap against, with the narrowPeak columns.

    Returns
    -------
    pandas.DataFrame
        Distinct rows of ``df1`` that pass.
    """
    if df1.empty or df2.empty:
        return df1.iloc[:0]

    overlaps = bf.overlap(df1, df2, return_overlap=True, how="inner")

    if overlaps.empty:
        return df1.iloc[:0]

    overlaps["len1"] = overlaps["end"] - overlaps["start"]
    overlaps["len2"] = overlaps["end_"] - overlaps["start_"]
    overlaps["overlap_len"] = overlaps["overlap_end"] - overlaps["overlap_start"]
    overlaps["frac1"] = overlaps["overlap_len"] / overlaps["len1"]
    overlaps["frac2"] = overlaps["overlap_len"] / overlaps["len2"]

    mask = (overlaps["frac1"] >= 0.5) | (overlaps["frac2"] >= 0.5)
    filtered = overlaps[mask]

    return filtered[PEAK_COLS].drop_duplicates()


def main(argv=None):
    """Filter the full-sample peaks against each pseudoreplicate in turn.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.

    Raises
    ------
    ValueError
        If no peak passes the first or the second pseudoreplicate.
    """
    args = build_parser().parse_args(argv)

    peaks = pd.read_csv(args.peaks, sep="\t", names=PEAK_COLS, header=None)
    pseudorep1 = pd.read_csv(args.pseudorep1, sep="\t", names=PEAK_COLS, header=None)
    pseudorep2 = pd.read_csv(args.pseudorep2, sep="\t", names=PEAK_COLS, header=None)

    passed_peaks = filter_50percent_overlap(peaks, pseudorep1)

    if passed_peaks.empty:
        raise ValueError("No peaks that overlap subset 1 by ≥50%")

    final_peaks = filter_50percent_overlap(passed_peaks, pseudorep2)

    if final_peaks.empty:
        raise ValueError("No peaks that overlap both subsets by ≥50%")

    final_peaks = final_peaks.drop_duplicates().sort_values(["chrom", "start"])
    final_peaks.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"Found {len(final_peaks)} peaks that overlap both subsets by ≥50%")


if __name__ == "__main__":
    main()
