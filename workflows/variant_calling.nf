include { EXTRACT_REGION } from '../modules/extract_region.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'
include { SORT_MARKDUP } from '../modules/sort_markdup.nf'
include { GATK_CALL } from '../modules/gatk_call.nf'
include { DEEPVARIANT_CALL } from '../modules/deepvariant_call.nf'

workflow VARIANT_CALLING {
    // Value channels (not Channel.fromPath queue channels) for anything
    // consumed by more than one process downstream — a queue channel drains
    // on its first use, so regions_bed/reference_fasta/reference_fai (each
    // needed by both an earlier and a later module) would come up empty the
    // second time otherwise.
    regions_bed = Channel.value(file(params.regions_bed))

    EXTRACT_REGION(regions_bed, params.hg002_bam_url)

    reference_fasta = Channel.value(file(params.reference_fasta))
    reference_fai    = Channel.value(file("${params.reference_fasta}.fai"))
    reference_dict   = Channel.value(file(params.reference_dict))
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

    SORT_MARKDUP(BWA_ALIGN.out.sam)

    GATK_CALL(
        SORT_MARKDUP.out.bam,
        SORT_MARKDUP.out.bai,
        reference_fasta,
        reference_fai,
        reference_dict,
        regions_bed
    )

    DEEPVARIANT_CALL(
        SORT_MARKDUP.out.bam,
        SORT_MARKDUP.out.bai,
        reference_fasta,
        reference_fai,
        regions_bed
    )

    // Next modules (batch 6+): cross_check_vcfs, happy_benchmark — see
    // project5_scoping.md pipeline architecture diagram.
}
