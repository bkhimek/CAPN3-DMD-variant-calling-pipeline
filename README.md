# CAPN3/DMD Variant Calling Pipeline

FASTQ → aligned BAM → called VCF, scoped to the CAPN3 and DMD loci, benchmarked
against a GIAB HG002 truth set. Upstream companion to
[CAPN3-DMD-variant-classifier](https://github.com/bkhimek/CAPN3-DMD-variant-classifier) —
this pipeline's VCF output is the intended eventual feed into that project's
`VariantEvidenceBundle` input.

See `project5_scoping.md` for the full design rationale and decisions, and
`docs/` for pinned data sources and coordinate derivation.

## Status

Batch 1 of N. `EXTRACT_REGION` (module 1 of 7) is implemented and verified
end-to-end: pulls the CAPN3 + DMD padded regions directly out of the remote
GIAB HG002 60x GRCh38 BAM via HTTP range requests (no full-genome download),
converts to paired FASTQ. Confirmed 336,549 paired reads extracted across
both loci.

Remaining modules (not yet built): `bwa_align`, `sort_markdup`, `gatk_call`,
`deepvariant_call`, `cross_check_vcfs`, `happy_benchmark`.

## Pipeline architecture

```
HG002 public pre-aligned BAM (GRCh38, remote, region-extracted via HTTP range + .bai)
        ▼
  region-subset BAM → FASTQ (R1/R2)          [EXTRACT_REGION — done]
        │  BWA-MEM2 align to GRCh38
        ▼
  aligned BAM → sort + mark duplicates
        ├──────────────┐
        ▼              ▼
  GATK HaplotypeCaller  DeepVariant
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

```bash
conda activate nextflow
nextflow run main.nf -profile docker
```

Outputs land in `results/` (gitignored).

## Repo conventions

Follows this workspace's Track A Nextflow DSL2 conventions (see
`~/projects/CLAUDE.md`): one process per module file in `modules/`, workflow
composition in `workflows/`, Docker-only tooling (no conda/Nextflow-managed
tool installs), resource limits capped at 4 CPU / 8GB per process.
