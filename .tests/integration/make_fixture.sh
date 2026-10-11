#!/usr/bin/env bash
# Builds resources/fixture: chr21 and chrM reads of two ENCODE ATAC-seq replicates,
# the chr21 + chrM genome, and the ENCODE references cut to those contigs.
# S1 and S10 are ENCSR095QNB replicates 1 and 2. Their IDs prefix one another on purpose.
# Needs the Weng Lab ENCODE mirror. Run from the repository root:
#   pixi run --manifest-path workflow/envs/pixi.toml -e test bash .tests/integration/make_fixture.sh
set -euo pipefail

out=resources/fixture
encode=/data/projects/encode/data/ENCSR095QNB
refs=/data/zusers/ramirezc/atac-smk-refs/GRCh38
fasta=/data/zusers/ramirezc/atac-smk/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
mkdir -p "$out"

samtools faidx "$fasta" chr21 chrM > "$out/genome.fa"

for accession in ENCFF493CCB ENCFF356LFX ENCFF304XEX ENCFF140XLU ENCFF212UAV; do
    zcat "$refs/$accession.bed.gz" \
    | awk '$1 == "chr21" || $1 == "chrM"' \
    | gzip -n \
    > "$out/$accession.bed.gz"
done

# Keeps pairs whose first mate ENCODE placed on chr21, and every tenth on chrM.
while read -r sample bam r1 r2; do
    samtools view -@ 8 -F 0x900 -f 64 "$encode/$bam.bam" \
    | awk '$3 == "chr21" || ($3 == "chrM" && mito++ % 10 == 0) { print $1 }' \
    > "$out/$sample.names"
    seqkit grep -j 8 -f "$out/$sample.names" "$encode/$r1.fastq.gz" -o "$out/${sample}_R1.fastq.gz"
    seqkit grep -j 8 -f "$out/$sample.names" "$encode/$r2.fastq.gz" -o "$out/${sample}_R2.fastq.gz"
    rm "$out/$sample.names"
done <<'EOF'
S1 ENCFF290DWH ENCFF411MVZ ENCFF736ZWY
S10 ENCFF878MKE ENCFF196VAZ ENCFF195FRB
EOF
