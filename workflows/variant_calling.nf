include { EXTRACT_REGION } from '../modules/extract_region.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'

workflow VARIANT_CALLING {
    regions_bed = Channel.fromPath(params.regions_bed)

    EXTRACT_REGION(regions_bed, params.hg002_bam_url)

    reference_fasta = Channel.fromPath(params.reference_fasta)
    reference_fai    = Channel.fromPath("${params.reference_fasta}.fai")
    reference_0123   = Channel.fromPath("${params.reference_fasta}.0123")
    reference_amb    = Channel.fromPath("${params.reference_fasta}.amb")
    reference_ann    = Channel.fromPath("${params.reference_fasta}.ann")
    reference_bwt    = Channel.fromPath("${params.reference_fasta}.bwt.2bit.64")
    reference_pac    = Channel.fromPath("${params.reference_fasta}.pac")

    BWA_ALIGN(
        EXTRACT_REGION.out.r1,
        EXTRACT_REGION.out.r2,
        reference_fasta,
        reference_fai,
        reference_0123,
        reference_amb,
        reference_ann,
        reference_bwt,
        reference_pac
    )

    // Next modules (batch 3+): sort_markdup, gatk_call, deepvariant_call,
    // cross_check_vcfs, happy_benchmark — see project5_scoping.md pipeline
    // architecture diagram.
}
