# Target region derivation

Verified live against Ensembl REST and NCBI Gene on 2026-07-28 (not taken from
the scoping doc's table, which flagged its own numbers as unverified).

## CAPN3 (ENSG00000092529 / NCBI Gene 825)

| Source | Assembly | Coordinates (1-based) |
|--------|----------|------------------------|
| Ensembl REST (`rest.ensembl.org/lookup/symbol`) | GRCh38 | chr15:42,359,494-42,412,949 |
| NCBI Gene (`NC_000015.10`) | GRCh38 | chr15:42,359,501-42,412,317 |

Union of both: chr15:42,359,494-42,412,949.

## DMD (ENSG00000198947 / NCBI Gene 1756)

| Source | Assembly | Coordinates (1-based) |
|--------|----------|------------------------|
| Ensembl REST | GRCh38 | chrX:31,097,677-33,339,609 |
| NCBI Gene (`NC_000023.11`) | GRCh38 | chrX:31,119,222-33,339,388 |

Union of both: chrX:31,097,677-33,339,609.

## Padding

+/-100kb per the scoping doc's rationale (mate-pair/split-read context near
gene boundaries). Final regions, 1-based inclusive:

- CAPN3: chr15:42,259,494-42,512,949
- DMD: chrX:30,997,677-33,439,609

`data/reference/regions.bed` encodes these in 0-based BED convention
(start - 1).

## Reference build note

The HG002 source BAM (see `docs/data_sources.md`) is aligned to
`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set` with UCSC-style
`chr`-prefixed contigs. All region coordinates above use that same
`chr`-prefixed convention (not Ensembl's bare `15`/`X` naming) so they line up
directly with the alignment reference and the GIAB truth VCF without a
liftover/rename step.
