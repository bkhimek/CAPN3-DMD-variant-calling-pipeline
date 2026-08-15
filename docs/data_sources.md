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

## Docker image tags (pinned 2026-07-28 unless noted per-row; verify again before a future rebuild)

| Tool | Image:tag |
|------|-----------|
| BWA-MEM2 | `quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0` |
| samtools | `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` |
| bcftools | `quay.io/biocontainers/bcftools:1.24--h487d631_1` |
| GATK | `broadinstitute/gatk:4.6.2.0` (also used for `SORT_MARKDUP` — confirmed this image bundles `samtools` at `/usr/bin/samtools` alongside `gatk`, so sort + `MarkDuplicates` + index share one container instead of pulling in a separate Picard-only image) |
| DeepVariant | `google/deepvariant:1.10.0` |
| hap.py | `jmcdani20/hap.py:v0.3.12` (pkrusche/hap.py is dead since 2017; this fork is the maintained alternative the scoping doc called out) |
| VEP | `ensemblorg/ensembl-vep:release_116.0` (pinned 2026-07-29; verified live against Docker Hub's tag list — `release_117.0` does not yet exist, `release_116.0` and `latest` share the same digest, both last pushed 2026-06-10) |
| FastQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` (pinned 2026-08-13; quay.io itself was unreachable from the sandbox this tag was checked from, so verified indirectly via nf-core's own current module source, which quotes this exact tag — FastQC hasn't had a release since 2023, consistent with this being current) |
| MultiQC | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1` (pinned 2026-08-13; live-verified against the quay.io tag API directly — 1.35 is MultiQC's current release, and `_1` is the newer of two builds published for it) |
| fastp | `quay.io/biocontainers/fastp:1.3.6--h43da1c4_0` (pinned 2026-08-14; live-verified against the quay.io tag API directly, matching fastp's current release, 1.3.6. Worth a note: two indirect sources checked first — a registry mirror and nf-core's own current module — both showed fastp's biocontainers image capping out a full major version earlier at `1.1.0`, and a `1.3.6--h43da1c4_0` tag quoted in a third-party blog post was initially rejected as probably pattern-guessed. The live quay.io API proved both indirect sources were simply stale, and the blog's tag was correct all along — the tag existed, it just hadn't propagated to either secondary source yet.) |

## Ensembl GRCh38 gene annotation (GFF3) — for `ANNOTATE_CALLS`/VEP

```
https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz
BSD sum (Ensembl's own CHECKSUMS file format): 02333 105838
```

Plain gzip, not bgzip/tabix-indexed at the source (confirmed live: no `.tbi`
alongside it), so — same rationale as the reference FASTA above — this is a
one-time ~100MB download via `scripts/fetch_gene_annotation.sh`, not fetched
per pipeline run. The script extracts the full CAPN3 (`ENSG00000092529`) and
DMD (`ENSG00000198947`) gene hierarchies (gene → transcript → exon/CDS/UTR)
by following GFF3 Parent references three levels deep, *not* a coordinate-
window slice: an earlier attempt that kept any line overlapping the padded
regions left dangling child features whose parent gene fell just outside the
window, which crashed VEP (`Can't call method "strand" on an undefined
value`) — caught and fixed before finalizing the script, not shipped broken.
Chromosome names are renamed from Ensembl's bare `15`/`X` to `chr15`/`chrX`
to match this project's BWA-aligned reference and VCFs. Output:
`data/reference/CAPN3_DMD.gff3.gz` (+ `.tbi`), ~40KB, gitignored.

MANE Select transcripts used for VEP's clinically-relevant-transcript
selection (cross-checked against NCBI/EBI's MANE summary file directly, not
carried over from memory): CAPN3 `NM_000070.3` / `ENST00000397163.8`; DMD
`NM_004006.3` / `ENST00000357033.9`.

## gnomAD v4.1 genomes — population frequency for `ANNOTATE_CALLS`

```
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr15.vcf.bgz(.tbi)
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrX.vcf.bgz(.tbi)
```

Confirmed live: real bgzip (`Accept-Ranges: bytes`, has a `.tbi`), so the same
remote-region-extraction trick used for the HG002 BAM and GIAB truth VCF
applies directly — `GNOMAD_ANNOTATE` pulls just the padded CAPN3/DMD regions
(one query per chromosome) from the 16.9GB chr15 file and the similarly large
chrX file, no full-genome download. `AF_grpmax` (max allele frequency across
genetic ancestry groups — gnomAD v4's field name for what used to be called
"popmax AF") is the source for `PopulationEvidence.ancestry_specific_max_af`;
confirmed via the source VCF's own INFO header, not assumed from an older
gnomAD version's field-naming convention.
