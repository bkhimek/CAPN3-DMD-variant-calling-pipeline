# CAPN3/DMD Variant Calling Pipeline

FASTQ → aligned BAM → called VCF, scoped to the CAPN3 and DMD loci, benchmarked
against a GIAB HG002 truth set. Upstream companion to
[CAPN3-DMD-variant-classifier](https://github.com/bkhimek/CAPN3-DMD-variant-classifier) —
this pipeline's VCF output is the intended eventual feed into that project's
`VariantEvidenceBundle` input.

See `project5_scoping.md` for the full design rationale and decisions, and
`docs/` for pinned data sources and coordinate derivation.

## Status

**Complete — core 7-module FASTQ→VCF→benchmark pipeline, plus an
`ANNOTATE_CALLS` annotation stage (2 modules) feeding
[CAPN3-DMD-variant-classifier](https://github.com/bkhimek/CAPN3-DMD-variant-classifier)'s
`VariantEvidenceBundle`.** All 9 modules implemented and verified
end-to-end:

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
- `CROSS_CHECK_VCFS` — `bcftools isec` concordance between `GATK_CALL` and
  `DEEPVARIANT_CALL`: 2,271 concordant calls, 47 GATK-only, and — after
  digging past the raw isec output — only 94 genuine DeepVariant-only
  variant calls. The raw isec count for DeepVariant-only records is 899, but
  805 of those are `RefCall`/`NoCall` sites: DeepVariant emits a record for
  every candidate site it examines, not only confident variant calls, so
  most of the raw private-record count isn't real caller disagreement. The
  module's `concordance_summary.txt` reports both numbers so this isn't
  silently mischaracterized. (Sanity check: 2,318 = 47 + 2,271 and 3,170 =
  899 + 2,271, both exact.)
- `HAPPY_BENCHMARK` — benchmarks each caller's VCF against the GIAB HG002
  NISTv4.2.1 truth set (remote region-extracted the same way as the source
  BAM, via `FETCH_TRUTH_SET`), restricted to the padded CAPN3/DMD regions,
  `--pass-only`. Both callers hit **perfect recall and precision (1.0/1.0)**
  for SNPs and INDELs in this region: GATK 545/545 SNP + 67/67 INDEL truth
  variants recovered with 0 false positives; DeepVariant 545/545 SNP +
  67/67 INDEL, also 0 false positives (DeepVariant's slightly higher
  QUERY.TOTAL reflects its extra non-PASS RefCall/NoCall records, correctly
  excluded by `--pass-only`, consistent with the `CROSS_CHECK_VCFS`
  finding). Both benchmark runs report the same `TRUTH.TOTAL` (612),
  confirming both correctly used the same fetched truth-set region.

- `VEP_ANNOTATE` — runs Ensembl VEP against the concordant (GATK ∩
  DeepVariant) call set — this pipeline's highest-confidence answer — using a
  gene-scoped GFF3 (CAPN3 + DMD full gene hierarchy, see
  `scripts/fetch_gene_annotation.sh`) and the existing scoped reference FASTA,
  rather than a full VEP cache download. Confirmed live: CSQ annotations on
  both genes' MANE Select transcripts (CAPN3 `ENST00000397163.8`, DMD
  `ENST00000357033.9`) are present and internally consistent (e.g. a CAPN3
  `c.706G>A` missense call's codon change `Gca/Aca` correctly predicts
  `p.Ala236Thr`; a DMD `c.8810G>A` call on the minus strand shows genomic
  `C>T`, consistent with reverse-complementing across the gene's strand).
- `GNOMAD_ANNOTATE` — annotates VEP's output with gnomAD v4.1 genomes
  population frequency (AC/AN/AF/AF_grpmax/nhomalt), region-extracted
  directly from the remote public sites VCFs (same HTTP-range + tabix trick
  as the HG002 BAM and GIAB truth VCF elsewhere in this pipeline — no
  16.9GB/chromosome full download). Confirmed live against an independent
  `bcftools view` lookup on two known concordant variants: CAPN3
  `chr15:42389001 G>A` (gnomAD AF=0.222151, matches exactly) and DMD
  `chrX:31478233 C>T` (gnomAD AF=0.938899, matches exactly) — both common
  variants, consistent with this pipeline's earlier finding that HG002's
  CAPN3/DMD calls are ordinary benign background polymorphism, not disease
  variants (see `~/projects/HANDOFF.md`).

All 9 modules done — the pipeline runs FASTQ (remote) → BAM → VCF →
benchmarked precision/recall report → annotated (transcript consequence +
population frequency) call set end-to-end.

Before running, fetch the scoped reference once (see `scripts/fetch_reference.sh`
below) and the gene annotation once (see `scripts/fetch_gene_annotation.sh`)
— both gitignored (large/regeneratable), not committed.

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
       cross-check (concordant / discordant calls)  [CROSS_CHECK_VCFS — done]
               ├──────────────────────────────┐
               ▼                              ▼
   hap.py vs. GIAB HG002 truth VCF      VEP (transcript consequence)
   [HAPPY_BENCHMARK — done]             [VEP_ANNOTATE — done]
               ▼                              ▼
     precision/recall report      gnomAD v4.1 (population frequency)
     per caller, per region       [GNOMAD_ANNOTATE — done]
                                          ▼
                              annotated_calls.vcf.gz
                       (→ CAPN3-DMD-variant-classifier's VariantEvidenceBundle adapter)
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

...and fetch the gene-scoped annotation for VEP (~40KB on disk once built,
verified against Ensembl's published checksum):

```bash
./scripts/fetch_gene_annotation.sh
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
