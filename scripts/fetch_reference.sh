#!/usr/bin/env bash
# One-time setup: fetches the exact GRCh38 reference the HG002 source BAM was
# aligned against, verifies it, and extracts chr15+chrX (full chromosomes —
# not just the padded CAPN3/DMD region — so alignment coordinates stay in
# true genome space, comparable to the GIAB truth VCF). Not run as part of
# the Nextflow pipeline: an 875MB download decompressing to ~3.1GB is wasteful
# to repeat on every pipeline invocation, and Nextflow's work/ cache is
# ephemeral. Run this once; the pipeline then reads the scoped output.
#
# Source: docs/data_sources.md
set -euo pipefail

REF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/reference"
CACHE_DIR="${REF_DIR}/.cache"
mkdir -p "$CACHE_DIR"

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
EXPECTED_MD5="a056c57649f3c9964c68aead3849bbf8"
GZ_PATH="${CACHE_DIR}/GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
FNA_PATH="${CACHE_DIR}/GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
OUT_FASTA="${REF_DIR}/GRCh38_chr15_chrX.fna"

if [[ ! -f "$GZ_PATH" ]]; then
    echo "Downloading reference (~875MB)..."
    curl -L --fail -o "$GZ_PATH" "$URL"
fi

echo "Verifying checksum..."
ACTUAL_MD5=$(md5sum "$GZ_PATH" | awk '{print $1}')
if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "MD5 mismatch: expected $EXPECTED_MD5, got $ACTUAL_MD5" >&2
    exit 1
fi
echo "Checksum OK ($ACTUAL_MD5)."

if [[ ! -f "$FNA_PATH" ]]; then
    echo "Decompressing..."
    gunzip -k -c "$GZ_PATH" > "$FNA_PATH"
fi

echo "Indexing full reference..."
docker run --rm -v "${CACHE_DIR}:/data" quay.io/biocontainers/samtools:1.24--h9dcdb79_1 \
    samtools faidx "/data/$(basename "$FNA_PATH")"

echo "Extracting chr15 + chrX..."
docker run --rm -v "${CACHE_DIR}:/data" -v "${REF_DIR}:/out" quay.io/biocontainers/samtools:1.24--h9dcdb79_1 \
    samtools faidx "/data/$(basename "$FNA_PATH")" chr15 chrX -o "/out/$(basename "$OUT_FASTA")"

echo "Indexing scoped reference..."
docker run --rm -v "${REF_DIR}:/out" quay.io/biocontainers/samtools:1.24--h9dcdb79_1 \
    samtools faidx "/out/$(basename "$OUT_FASTA")"

echo "Building BWA-MEM2 index..."
docker run --rm -v "${REF_DIR}:/out" quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0 \
    bwa-mem2 index "/out/$(basename "$OUT_FASTA")"

echo "Building GATK sequence dictionary..."
DICT_PATH="${REF_DIR}/$(basename "${OUT_FASTA%.fna}").dict"
rm -f "$DICT_PATH"
docker run --rm -v "${REF_DIR}:/out" broadinstitute/gatk:4.6.2.0 \
    gatk CreateSequenceDictionary -R "/out/$(basename "$OUT_FASTA")" -O "/out/$(basename "$DICT_PATH")"

echo "Done. Scoped reference at: $OUT_FASTA"
echo "Full-genome cache kept at: $CACHE_DIR (safe to delete to reclaim ~4GB; re-run this script to regenerate)"
