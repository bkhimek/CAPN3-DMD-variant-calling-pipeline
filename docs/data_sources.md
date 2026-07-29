# Pinned data sources

Verified reachable (HTTP 200, byte-range support checked where relevant) on
2026-07-28.

## HG002 pre-aligned GRCh38 BAM

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam
(+ .bam.bai alongside)
```

- ~60x coverage (first 20% of the 300x Illumina HiSeq 2x148bp AJ trio run),
  novoalign v3.02.07 against `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set`
  (UCSC-style `chr`-prefixed contigs).
- Chose the 60x file over the 300x file (also present in the same directory)
  — standard germline SNV/indel coverage, and both support HTTP range
  requests (`Accept-Ranges: bytes` confirmed) so `samtools view <url> <region>`
  can pull just the CAPN3/DMD loci via the `.bai` index without downloading
  the full ~126GB (300x) / correspondingly smaller-but-still-large 60x file.

## Reference genome — GRCh38, chr15 + chrX

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
MD5: a056c57649f3c9964c68aead3849bbf8 (from the directory's md5checksums.txt)
```

The exact reference the HG002 source BAM was aligned against (per its own
`README_NHGRI_Novoalign_bams`) — required so alignment coordinates and
variant calls come out in the same coordinate system as both the source BAM
and the GIAB truth VCF.

This file is plain gzip, not bgzip (confirmed: `samtools faidx` on the remote
URL fails with "Cannot index files compressed with gzip, please use bgzip"),
so unlike the HG002 BAM there's no free remote random-access trick here — it
has to be downloaded in full (~875MB compressed) once. `scripts/fetch_reference.sh`
downloads it, verifies the checksum above, decompresses, then extracts full
chr15 + chrX (not just the padded CAPN3/DMD region — alignment needs the full
chromosome as the reference contig so BWA-MEM2's reported positions are true
GRCh38 genome coordinates, not offsets from an arbitrary region start) and
builds the samtools + BWA-MEM2 indexes. Output: `data/reference/GRCh38_chr15_chrX.fna`
(+ indexes), ~1.6GB on disk, gitignored. The full-genome intermediate is
deleted after extraction to reclaim ~4GB; re-run the script to regenerate.

Verified: extracted chr15 is 101,991,189bp and chrX is 156,040,895bp — both
match the known GRCh38 chromosome lengths exactly.

## GIAB truth set — NISTv4.2.1, GRCh38

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz(.tbi)
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

Chose v4.2.1 over the newer v5.0q: v4.2.1 is a single small-variant
benchmark VCF+BED (matches this project's SNV/indel-only scope); v5.0q
splits into separate `smvar`/`stvar` benchmark sets and is CHM13-primary
with GRCh38 as a lifted-over secondary release — more moving parts than this
project needs. Revisit if project 5 later grows a structural-variant track.

## Docker image tags (pinned 2026-07-28, verify again before a future rebuild)

| Tool | Image:tag |
|------|-----------|
| BWA-MEM2 | `quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0` |
| samtools | `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` |
| bcftools | `quay.io/biocontainers/bcftools:1.24--h487d631_1` |
| GATK | `broadinstitute/gatk:4.6.2.0` (also used for `SORT_MARKDUP` — confirmed this image bundles `samtools` at `/usr/bin/samtools` alongside `gatk`, so sort + `MarkDuplicates` + index share one container instead of pulling in a separate Picard-only image) |
| DeepVariant | `google/deepvariant:1.10.0` |
| hap.py | `jmcdani20/hap.py:v0.3.12` (pkrusche/hap.py is dead since 2017; this fork is the maintained alternative the scoping doc called out) |
