"""Convert a filtered BAM to BEDPE, then to a sorted, gzipped tagAlign."""

import argparse
import os
import subprocess
import tempfile


def build_parser():
    """Build the bam_to_tagalign.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the bam_to_tagalign.py flags.
    """
    p = argparse.ArgumentParser(
        description="Convert BAM to BEDPE and then to sorted TagAlign format for ATAC-seq."
    )
    p.add_argument("--bam", required=True, help="Filtered BAM.")
    p.add_argument("--bedpe", required=True, help="BEDPE to write.")
    p.add_argument(
        "--tagalign",
        required=True,
        help="TagAlign to write; gzip replaces it with <tagalign>.gz.",
    )
    p.add_argument("--threads", type=int, default=1, help="Threads for samtools sort.")
    return p


def check_bam_sort_order(bam_file):
    """Check whether a BAM header declares queryname sort order.

    Parameters
    ----------
    bam_file : str
        BAM to check.

    Returns
    -------
    bool
        True if the ``@HD`` line has ``SO:queryname``.
    """
    cmd = f"samtools view -H {bam_file} | grep '@HD' | grep 'SO:queryname'"
    result = subprocess.run(cmd, shell=True, capture_output=True, check=False)
    return result.returncode == 0


def sort_bam_by_name(bam_file, sorted_bam, threads):
    """Sort a BAM by read name with samtools sort.

    Parameters
    ----------
    bam_file : str
        BAM to sort.
    sorted_bam : str
        Name-sorted BAM to write.
    threads : int
        Threads for samtools sort.
    """
    try:
        cmd = f"samtools sort -n -@ {threads} -o {sorted_bam} {bam_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Sorted BAM by name: {sorted_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error sorting BAM by name: {e}")


def bam_to_bedpe(bam_file, output_bedpe):
    """Convert a name-sorted BAM to BEDPE with bedtools.

    Parameters
    ----------
    bam_file : str
        Name-sorted BAM.
    output_bedpe : str
        BEDPE to write.
    """
    try:
        cmd = f"bedtools bamtobed -bedpe -i {bam_file} > {output_bedpe}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Converted BAM to BEDPE: {output_bedpe}")
    except subprocess.CalledProcessError as e:
        print(f"Error during BAM to BEDPE conversion: {e}")


def bedpe_to_tagalign(bedpe_file, output_tagalign):
    """Write each BEDPE mate as a Tn5-shifted tagAlign line, sorted by position.

    Parameters
    ----------
    bedpe_file : str
        BEDPE from ``bam_to_bedpe``.
    output_tagalign : str
        Sorted tagAlign to write.
    """
    try:
        with (
            tempfile.NamedTemporaryFile(
                "w", delete=False, suffix=".tagAlign"
            ) as tagalign,
            open(bedpe_file) as bedpe,
        ):
            for line in bedpe:
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    continue

                chrom1, start1, end1, chrom2, start2, end2, _, _, strand1, strand2 = (
                    fields
                )

                # Moves each read end to the Tn5 insertion site.
                if strand1 == "+":
                    start1 = int(start1) + 4
                elif strand1 == "-":
                    end1 = int(end1) - 5

                if strand2 == "+":
                    start2 = int(start2) + 4
                elif strand2 == "-":
                    end2 = int(end2) - 5

                tagalign.write(f"{chrom1}\t{start1}\t{end1}\tN\t1000\t{strand1}\n")
                tagalign.write(f"{chrom2}\t{start2}\t{end2}\tN\t1000\t{strand2}\n")

        print(f"Converted BEDPE to temporary TagAlign: {tagalign.name}")

        cmd = f"sort -k1,1 -k2,2n {tagalign.name} > {output_tagalign}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Sorted TagAlign file: {output_tagalign}")

        os.unlink(tagalign.name)

    except Exception as e:  # noqa: BLE001
        print(f"Error during BEDPE to TagAlign conversion: {e}")


def main(argv=None):
    """Name-sort the BAM if needed, then write the BEDPE and gzipped tagAlign.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    temp_bam = None
    bam_to_use = args.bam
    if not check_bam_sort_order(args.bam):
        print("BAM file is not sorted by name. Sorting...")
        fd, temp_bam = tempfile.mkstemp(suffix=".bam")
        os.close(fd)
        sort_bam_by_name(args.bam, temp_bam, args.threads)
        bam_to_use = temp_bam

    bam_to_bedpe(bam_to_use, args.bedpe)

    bedpe_to_tagalign(args.bedpe, args.tagalign)

    try:
        cmd = f"gzip -f {args.tagalign}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Compressed TagAlign file: {args.tagalign}.gz")
    except subprocess.CalledProcessError as e:
        print(f"Error during gzip compression: {e}")

    if temp_bam:
        os.unlink(temp_bam)
        print(f"Deleted temporary BAM file: {temp_bam}")


if __name__ == "__main__":
    main()
