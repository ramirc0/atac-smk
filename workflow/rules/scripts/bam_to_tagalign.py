import subprocess
import argparse
import os
import tempfile

def check_bam_sort_order(bam_file):
    """
    Check if a BAM file is sorted by name.
    
    :param bam_file: Path to the BAM file.
    :return: True if sorted by name, False otherwise.
    """
    try:
        cmd = f"samtools view -H {bam_file} | grep '@HD' | grep 'SO:queryname'"
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.returncode == 0  # Return True if BAM is sorted by name
    except subprocess.CalledProcessError as e:
        print(f"Error checking BAM sort order: {e}")
        return False

def sort_bam_by_name(bam_file, sorted_bam, threads):
    """
    Sort a BAM file by name using samtools sort.
    
    :param bam_file: Path to the input BAM file.
    :param sorted_bam: Path to the output sorted BAM file.
    :param threads: Number of threads to use for sorting.
    """
    try:
        cmd = f"samtools sort -n -@ {threads} -o {sorted_bam} {bam_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Sorted BAM by name: {sorted_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error sorting BAM by name: {e}")

def bam_to_bedpe(bam_file, output_bedpe):
    """
    Convert a BAM file to BEDPE format using bedtools.
    
    :param bam_file: Path to the input BAM file.
    :param output_bedpe: Path to the output BEDPE file.
    """
    try:
        cmd = f"bedtools bamtobed -bedpe -i {bam_file} > {output_bedpe}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Converted BAM to BEDPE: {output_bedpe}")
    except subprocess.CalledProcessError as e:
        print(f"Error during BAM to BEDPE conversion: {e}")

def bedpe_to_tagalign(bedpe_file, output_tagalign):
    """
    Convert a BEDPE file to TagAlign format with ATAC-seq Tn5 shift.
    
    :param bedpe_file: Path to the input BEDPE file.
    :param output_tagalign: Path to the output TagAlign file.
    """
    try:
        temp_tagalign = tempfile.NamedTemporaryFile(delete=False, suffix=".tagAlign")
        
        with open(bedpe_file, 'r') as bedpe, open(temp_tagalign.name, 'w') as tagalign:
            for line in bedpe:
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue  # Skip malformed lines
                
                chrom1, start1, end1, chrom2, start2, end2, _, _, strand1, strand2 = fields
                
                # Tn5 shift: adjust start and end positions
                if strand1 == "+":
                    start1 = int(start1) + 4
                elif strand1 == "-":
                    end1 = int(end1) - 5
                
                if strand2 == "+":
                    start2 = int(start2) + 4
                elif strand2 == "-":
                    end2 = int(end2) - 5

                # Write to TagAlign format (2 lines per fragment)
                tagalign.write(f"{chrom1}\t{start1}\t{end1}\tN\t1000\t{strand1}\n")
                tagalign.write(f"{chrom2}\t{start2}\t{end2}\tN\t1000\t{strand2}\n")
        
        print(f"Converted BEDPE to temporary TagAlign: {temp_tagalign.name}")

        # Sort the TagAlign file
        cmd = f"sort -k1,1 -k2,2n {temp_tagalign.name} > {output_tagalign}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Sorted TagAlign file: {output_tagalign}")

        # Delete temporary file
        os.unlink(temp_tagalign.name)

    except Exception as e:
        print(f"Error during BEDPE to TagAlign conversion: {e}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert BAM to BEDPE and then to sorted TagAlign format for ATAC-seq.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("bedpe_file", help="Output BEDPE file")
    parser.add_argument("tagalign_file", help="Output sorted TagAlign file (gzipped)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for samtools sorting (default: 1)")

    args = parser.parse_args()

    # Step 1: Check if BAM is sorted by name; if not, sort it
    temp_bam = None
    bam_to_use = args.bam_file
    if not check_bam_sort_order(args.bam_file):
        print("BAM file is not sorted by name. Sorting...")
        temp_bam = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
        sort_bam_by_name(args.bam_file, temp_bam.name, args.threads)
        bam_to_use = temp_bam.name

    # Step 2: Convert BAM to BEDPE
    bam_to_bedpe(bam_to_use, args.bedpe_file)

    # Step 3: Convert BEDPE to TagAlign with Tn5 shift and sort
    bedpe_to_tagalign(args.bedpe_file, args.tagalign_file)

    # Step 4: Gzip the TagAlign file
    try:
        cmd = f"gzip -f {args.tagalign_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Compressed TagAlign file: {args.tagalign_file}.gz")
    except subprocess.CalledProcessError as e:
        print(f"Error during gzip compression: {e}")

    # Step 5: Clean up temporary BAM file
    if temp_bam:
        os.unlink(temp_bam.name)
        print(f"Deleted temporary BAM file: {temp_bam.name}")

if __name__ == "__main__":
    main()
