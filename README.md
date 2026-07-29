# CAPN3/DMD Variant Calling Pipeline

FASTQ → aligned BAM → called VCF, scoped to the CAPN3 and DMD loci, benchmarked
against a GIAB HG002 truth set. Upstream companion to
[CAPN3-DMD-variant-classifier](https://github.com/bkhimek/CAPN3-DMD-variant-classifier) —
this pipeline's VCF output is the intended eventual feed into that project's
`VariantEvidenceBundle` input.

See `project5_scoping.md` for the full design rationale and decisions, and
`docs/` for pinned data sources and coordinate derivation.

## Status

Batch 5 of N. Modules 1-5 of 7 implemented and verified end-to-end:

- `EXTRACT_REGION` — pulls the CAPN3 + DMD padded regions directly out of the
  remote GIAB HG002 60x GRCh38 BAM via HTTP range requests (no full-genome
  download), converts to paired FASTQ. Confirmed 336,549 paired reads
  extracted across both loci.
- `BWA_ALIGN` — aligns those reads against a scoped GRCh38 reference (full
  chr15 + chrX chromosomes, not just the padded region, so alignment
  coordinates stay in true genome space and line up with the GIAB truth VCF).
  Confirmed 99.76% mapped, 99.22% properly paired, and every mapped read
  lands on chr15 or chrX (no off-target contamination).
- `SORT_MARKDUP` — coordinate-sorts the aligned reads and flags PCR/optical
  duplicates via GATK `MarkDuplicates` (the `broadinstitute/gatk` image
  bundles samtools too, so sort + dedup + index share one container — no
  separate Picard image needed). Confirmed flagstat totals unchanged from
  `BWA_ALIGN` (dedup flags reads, doesn't remove them) and a 0.95% duplication
  rate (6,370 of 673,098 primary reads), consistent with expected library
  complexity for this coverage.
- `GATK_CALL` — runs GATK `HaplotypeCaller` restricted to the padded
  CAPN3/DMD regions (`-L data/reference/regions.bed`), single-sample HG002.
  Confirmed 2,318 variants called (702 on chr15, 1,616 on chrX), every
  variant's position falling inside the padded gene regions, with sane
  QUAL/DP/genotype fields. Needed the reference sequence dictionary
  (`.dict`), which `scripts/fetch_reference.sh` now also builds via GATK
  `CreateSequenceDictionary`.
- `DEEPVARIANT_CALL` — runs DeepVariant (WGS model) on the same
  `SORT_MARKDUP` BAM, region-restricted the same way as `GATK_CALL`, for the
  later cross-caller check. Confirmed 3,170 variants called (837 chr15,
  2,333 chrX) — a different, larger count than GATK's is expected (different
  caller/algorithm), and the two callers' first records agree exactly
  (`chr15:42259883 C>T`), a good early sanity signal ahead of the dedicated
  `CROSS_CHECK_VCFS` module.

Remaining modules (not yet built): `cross_check_vcfs`, `happy_benchmark`.

Before running, fetch the scoped reference once (see `scripts/fetch_reference.sh`
below) — it's gitignored (large binary), not committed.

## Pipeline architecture

```
HG002 public pre-aligned BAM (GRCh38, remote, region-extracted via HTTP range + .bai)
        ▼
  region-subset BAM → FASTQ (R1/R2)          [EXTRACT_REGION — done]
        │  BWA-MEM2 align to GRCh38          [BWA_ALIGN — done]
        ▼
  aligned BAM → sort + mark duplicates       [SORT_MARKDUP — done]
        ├──────────────┐
        ▼              ▼
  GATK HaplotypeCaller  DeepVariant
  [GATK_CALL — done]    [DEEPVARIANT_CALL — done]
        └──────┬───────┘
               ▼
       cross-check (concordant / discordant calls)
               ▼
   hap.py vs. GIAB HG002 truth VCF (region-subset confident BED)
               ▼
     precision/recall report per caller, per region
```

## Target regions (GRCh38, chr-prefixed, verified 2026-07-28 — see `docs/regions.md`)

| Gene | Region (1-based) |
|------|-------------------|
| CAPN3 | chr15:42,259,494-42,512,949 |
| DMD | chrX:30,997,677-33,439,609 |

## Pinned data sources

See `docs/data_sources.md` for the exact HG002 BAM URL, GIAB truth set
version (NISTv4.2.1/GRCh38), and Docker image tags — all reachability- and
(where relevant) random-access-verified, not assumed.

## Running

One-time setup — fetch and index the scoped reference (chr15 + chrX,
~1.6GB on disk once built, verified against NCBI's published checksum):

```bash
./scripts/fetch_reference.sh
```

Then run the pipeline:

```bash
conda activate nextflow
nextflow run main.nf -profile docker
```

Outputs land in `results/` (gitignored). Use `-resume` to reuse cached
results from prior modules when iterating on a new one.

## Repo conventions

Follows this workspace's Track A Nextflow DSL2 conventions (see
`~/projects/CLAUDE.md`): one process per module file in `modules/`, workflow
composition in `workflows/`, Docker-only tooling (no conda/Nextflow-managed
tool installs), resource limits capped at 4 CPU / 8GB per process.
