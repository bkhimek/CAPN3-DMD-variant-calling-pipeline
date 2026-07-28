include { EXTRACT_REGION } from '../modules/extract_region.nf'

workflow VARIANT_CALLING {
    regions_bed = Channel.fromPath(params.regions_bed)

    EXTRACT_REGION(regions_bed, params.hg002_bam_url)

    // Next modules (batch 2+): bwa_align, sort_markdup, gatk_call,
    // deepvariant_call, cross_check_vcfs, happy_benchmark — see
    // project5_scoping.md pipeline architecture diagram.
}
