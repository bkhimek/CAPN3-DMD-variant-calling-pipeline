# Project 5 Scoping — Sequencing Reads → Variant Calling Pipeline

Status: scoped, not yet built. This document is meant to be pasted (alongside `CLAUDE.md`) into a new Claude Code session on WSL to kick off the build.

## Goal

FASTQ → aligned BAM → called VCF, benchmarked against a GIAB truth set, with the eventual VCF output feeding project 4's (CAPN3/DMD classifier) `VariantEvidenceBundle` input.

## Decisions

### Environment: WSL via Claude Code (not Cowork)

Project 4's Cowork sandbox pattern (stdlib Python, PyYAML, zip-and-sync) exists because that project needed no real bioinformatics tools. This project needs an aligner, two variant callers, samtools/bcftools, and Docker — and has to move FASTQ/BAM files that don't fit a zip-sync workflow. That profile matches Track A's WSL2 + conda/Docker/Nextflow setup, not Track B's sandbox. Build this one in a WSL Claude Code session with direct file access, following Track A's Nextflow DSL2 conventions.

### Tooling: BWA-MEM2 + DeepVariant + GATK, cross-checked

- **Aligner:** BWA-MEM2 — still the standard short-read aligner used in GIAB benchmarking papers (minimap2/minibwa are faster but less established for this use case).
- **Callers, both, cross-checked:** GATK HaplotypeCaller (most mature best-practices documentation, traditional "gold standard") and DeepVariant (higher precision/recall against GIAB truth and lower Mendelian error rate in recent trio studies — Google's neural-net caller, GIAB-native since it's benchmarked in the PrecisionFDA Truth Challenges). Running both and cross-checking mirrors project 4's own pattern of implementing two independent systems (classic Table 5 vs. Bayesian) and cross-validating them against each other.
- **Benchmarking tool:** hap.py (Illumina) — the field-standard tool for scoring a VCF against a GIAB truth VCF + confident-region BED. `pkrusche/hap.py` is the reference Docker image; `jmcdani20/hap.py` bundles rtg-tools's `vcfeval` as an alternative comparison engine if useful.

### Scope: targeted regions (CAPN3 + DMD loci), not whole genome

A full HG002 WGS run is ~300x PCR-free coverage, tens of GB of FASTQ, hours of align+call compute — a slow loop for a portfolio project. Restricting to the CAPN3 and DMD loci (padded) keeps the pipeline fast to iterate on, keeps repo/data size sane, and ties directly into project 4's genes. Trade-off: doesn't exercise genome-wide edge cases (repetitive regions elsewhere, whole-genome QC) — acceptable for this project's purpose (pipeline correctness demo + direct feed into project 4).

**Gene coordinates (GRCh38) — confirm exact values against the actual reference GTF used in the build, don't hardcode these blindly:**

| Gene | Locus | Approx. span (Ensembl) | Approx. span (NCBI RefSeq) |
|------|-------|------------------------|------------------------------|
| CAPN3 | chr15q15.1 | chr15:42,359,498–42,412,949 | chr15:42,359,501–42,412,317 |
| DMD | chrXp21.2-p21.1 | chrX:31,097,677–33,339,609 | chrX:31,119,222–33,339,388 |

Ensembl/RefSeq annotations disagree slightly at the boundaries (isoform differences) — re-derive from the GTF/GFF actually used as the pipeline reference rather than trusting this table at build time.

Recommend padding each locus by ~100 kb on each side to safely capture split reads / mate pairs near boundaries without ballooning data size (SNV/indel calling doesn't need more than typical insert-size context beyond the gene body).

## Pipeline architecture

```
HG002 public pre-aligned BAM/CRAM (GRCh38)
        │  samtools view -b <region>  (chr15 + chrX loci, padded)
        ▼
  region-subset BAM  →  samtools fastq  →  region FASTQ (R1/R2)
        │  BWA-MEM2 align to GRCh38
        ▼
  aligned BAM  →  sort + mark duplicates
        ├──────────────┐
        ▼              ▼
  GATK HaplotypeCaller  DeepVariant
        │              │
        ▼              ▼
     GATK VCF      DeepVariant VCF
        │              │
        └──────┬───────┘
               ▼
       cross-check (concordant / discordant calls)
               │
               ▼
   hap.py vs. GIAB HG002 truth VCF (region-subset via confident BED)
               │
               ▼
     precision/recall report per caller, per region
```

Extracting reads from a public pre-aligned HG002 BAM/CRAM (rather than downloading full raw FASTQ and aligning genome-wide) is the key trick that keeps this tractable — GIAB publishes pre-aligned GRCh38 BAMs (e.g. under the NHGRI_Illumina300X collection, also mirrored on public AWS/GCP buckets). Re-deriving FASTQ from the region-subset BAM keeps the "start from raw reads" framing honest while avoiding a multi-hundred-GB download.

## Repo structure (Track A Nextflow DSL2 conventions)

```
<repo-name>/
├── main.nf
├── nextflow.config          # docker + singularity profiles, executor limits
├── modules/
│   ├── extract_region.nf    # samtools view region subset + fastq conversion
│   ├── bwa_align.nf
│   ├── sort_markdup.nf
│   ├── gatk_call.nf
│   ├── deepvariant_call.nf
│   ├── cross_check_vcfs.nf
│   └── happy_benchmark.nf
├── workflows/
│   └── variant_calling.nf
├── data/
│   ├── reference/            # GRCh38 reference FASTA + index (or region-subset)
│   └── truth/                # GIAB HG002 truth VCF + confident-region BED (subset)
├── results/                  # gitignored
├── docs/
│   └── design notes specific to project 5
└── .gitignore
```

`.gitignore`: exclude `work/`, `results/`, `.nextflow/`, plus raw sequencing data (`*.bam`, `*.cram`, `*.fastq*`, `*.fna`, `*.fasta` except small curated fixtures) — consistent with both tracks' "never commit large real sequencing files" rule.

## Docker images (verify current tags at build time — registries move)

| Tool | Image (verify tag on Docker Hub at build time) |
|------|--------------------------------------------------|
| BWA-MEM2 | `quay.io/biocontainers/bwa-mem2:<latest>` |
| samtools/bcftools | `quay.io/biocontainers/samtools:<latest>` |
| GATK | `broadinstitute/gatk:<latest stable>` (≈4.5–4.6.x as of mid-2026) |
| DeepVariant | `google/deepvariant:<latest>` (1.9.0 was the newest numbered release found) |
| hap.py | `pkrusche/hap.py` (or `jmcdani20/hap.py` for bundled rtg-tools) |

## Open items for the WSL kickoff session

1. Confirm exact CAPN3/DMD coordinates against the reference GTF actually used.
2. Locate and pin a specific public HG002 pre-aligned GRCh38 BAM/CRAM URL (GIAB FTP / AWS / GCP mirror) and a specific GIAB truth-set version (v4.2.1 is the stable widely-used one; v5.0q is newer if preferred).
3. Confirm current stable Docker tags for each tool above.
4. Decide repo name (suggest: `HG002-variant-calling-pipeline`) and create it under `~/projects/`.
5. Update this project's entry in `CLAUDE.md` and start `HANDOFF.md` once underway.

## Verification discipline

Same standard as project 4: every claim gets checked against a known-correct source, nothing is asserted from memory. Concretely — hap.py's precision/recall numbers against GIAB truth are the ground-truth check for the pipeline itself, and GATK/DeepVariant concordance is a secondary internal cross-check (agreement isn't proof of correctness, but disagreement flags calls worth inspecting by hand before trusting them).
