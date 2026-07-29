#!/usr/bin/env bash
# One-time setup: fetches Ensembl's full GRCh38 GFF3 gene annotation (plain
# gzip, not bgzip/tabix-indexed at the source — confirmed live, so unlike the
# HG002 BAM/GIAB truth VCF/gnomAD sites VCFs, there's no free remote-range
# trick here), then extracts just the CAPN3 and DMD gene structures (gene ->
# transcript -> exon/CDS/UTR, full hierarchy, not a blind coordinate slice —
# see the ANNOTATE_CALLS entry in docs/data_sources.md for why a naive
# coordinate-window filter over neighboring genes crashes VEP with dangling
# Parent references) into a small bgzip+tabix-indexed GFF3 subset for VEP's
# --gff mode. Not run as part of the Nextflow pipeline, same rationale as
# fetch_reference.sh: a ~100MB+ full-genome-annotation download is wasteful
# to repeat on every pipeline invocation.
#
# Source: docs/data_sources.md
set -euo pipefail

REF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/reference"
CACHE_DIR="${REF_DIR}/.cache"
mkdir -p "$CACHE_DIR"

RELEASE="116"
URL="https://ftp.ensembl.org/pub/release-${RELEASE}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${RELEASE}.gff3.gz"
EXPECTED_SUM="02333 105838"
GZ_PATH="${CACHE_DIR}/Homo_sapiens.GRCh38.${RELEASE}.gff3.gz"
OUT_GFF="${REF_DIR}/CAPN3_DMD.gff3.gz"

# CAPN3 = ENSG00000092529, DMD = ENSG00000198947 (confirmed against the MANE
# summary file cross-checked in docs/regions.md's companion annotation notes).
CAPN3_GENE="gene:ENSG00000092529"
DMD_GENE="gene:ENSG00000198947"

if [[ ! -f "$GZ_PATH" ]]; then
    echo "Downloading Ensembl release ${RELEASE} GFF3 (~100MB)..."
    curl -sL --fail -o "$GZ_PATH" "$URL"
fi

echo "Verifying checksum (BSD sum, per Ensembl's own CHECKSUMS file)..."
ACTUAL_SUM=$(sum "$GZ_PATH" | awk '{print $1, $2}')
if [[ "$ACTUAL_SUM" != "$EXPECTED_SUM" ]]; then
    echo "Checksum mismatch: expected '$EXPECTED_SUM', got '$ACTUAL_SUM'" >&2
    exit 1
fi
echo "Checksum OK ($ACTUAL_SUM)."

FULL_GFF="${CACHE_DIR}/Homo_sapiens.GRCh38.${RELEASE}.gff3"
if [[ ! -f "$FULL_GFF" ]]; then
    echo "Decompressing..."
    gunzip -k -c "$GZ_PATH" > "$FULL_GFF"
fi

WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

echo "Extracting CAPN3/DMD gene hierarchy (gene -> transcript -> exon/CDS/UTR)..."
# Step 1: the two gene features themselves.
awk -F'\t' -v g1="$CAPN3_GENE" -v g2="$DMD_GENE" \
    '!/^#/ && $3=="gene" && ($9 ~ ("ID="g1";") || $9 ~ ("ID="g2";"))' \
    "$FULL_GFF" > "${WORK_DIR}/step1_genes.gff3"

# Step 2: direct children of those genes (transcripts/mRNAs).
awk -F'\t' -v g1="$CAPN3_GENE" -v g2="$DMD_GENE" \
    '!/^#/ && ($9 ~ ("Parent="g1";") || $9 ~ ("Parent="g1"$") || $9 ~ ("Parent="g2";") || $9 ~ ("Parent="g2"$"))' \
    "$FULL_GFF" > "${WORK_DIR}/step2_children.gff3"

grep -oE 'ID=transcript:[A-Z0-9]+' "${WORK_DIR}/step2_children.gff3" | sed 's/ID=//' | sort -u > "${WORK_DIR}/transcript_ids.txt"

# Step 3: grandchildren of those genes (exon/CDS/UTR rows, keyed by transcript Parent).
awk -F'\t' 'NR==FNR{ids[$1]=1; next} !/^#/ { match($9,/Parent=[^;]+/); p=substr($9,RSTART+7,RLENGTH-7); if (p in ids) print }' \
    "${WORK_DIR}/transcript_ids.txt" "$FULL_GFF" > "${WORK_DIR}/step3_grandchildren.gff3"

{
    echo "##gff-version 3"
    cat "${WORK_DIR}/step1_genes.gff3" "${WORK_DIR}/step2_children.gff3" "${WORK_DIR}/step3_grandchildren.gff3" \
        | awk -F'\t' '{ $1 = ($1=="15" ? "chr15" : ($1=="X" ? "chrX" : $1)); print }' OFS='\t' \
        | sort -k1,1 -k4,4n
} > "${WORK_DIR}/CAPN3_DMD.gff3"

N_LINES=$(wc -l < "${WORK_DIR}/CAPN3_DMD.gff3")
echo "Extracted ${N_LINES} feature lines."

echo "Sorting/bgzipping/indexing..."
docker run --rm -v "${WORK_DIR}:/data" -w /data quay.io/biocontainers/bcftools:1.24--h487d631_1 \
    bash -c "bgzip -f CAPN3_DMD.gff3 && tabix -p gff CAPN3_DMD.gff3.gz"

cp "${WORK_DIR}/CAPN3_DMD.gff3.gz" "$OUT_GFF"
cp "${WORK_DIR}/CAPN3_DMD.gff3.gz.tbi" "${OUT_GFF}.tbi"

echo "Done. Gene-scoped annotation at: $OUT_GFF"
echo "Full-genome cache kept at: $CACHE_DIR (safe to delete to reclaim ~1GB decompressed; re-run this script to regenerate)"
