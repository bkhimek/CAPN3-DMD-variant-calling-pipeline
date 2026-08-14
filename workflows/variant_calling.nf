include { EXTRACT_REGION } from '../modules/extract_region.nf'
include { FASTQC } from '../modules/fastqc.nf'
include { MULTIQC } from '../modules/multiqc.nf'
include { TRIM_READS } from '../modules/trim_reads.nf'
include { BWA_ALIGN } from '../modules/bwa_align.nf'
include { SORT_MARKDUP } from '../modules/sort_markdup.nf'
include { GATK_CALL } from '../modules/gatk_call.nf'
include { DEEPVARIANT_CALL } from '../modules/deepvariant_call.nf'
include { CROSS_CHECK_VCFS } from '../modules/cross_check_vcfs.nf'
include { FETCH_TRUTH_SET } from '../modules/fetch_truth_set.nf'
include { HAPPY_BENCHMARK } from '../modules/happy_benchmark.nf'
include { VEP_ANNOTATE } from '../modules/vep_annotate.nf'
include { GNOMAD_ANNOTATE } from '../modules/gnomad_annotate.nf'

workflow VARIANT_CALLING {
    // Value channels (not Channel.fromPath queue channels) for anything
    // consumed by more than one process downstream — a queue channel drains
    // on its first use, so regions_bed/reference_fasta/reference_fai (each
    // needed by both an earlier and a later module) would come up empty the
    // second time otherwise.
    regions_bed = Channel.value(file(params.regions_bed))

    EXTRACT_REGION(regions_bed, params.hg002_bam_url)

    // Raw-read QC/trimming extension, batch 1: FastQC + MultiQC on the raw
    // EXTRACT_REGION output, before any trimming. TRIM_READS (batch 2) and a
    // second FASTQC/MULTIQC pass on trimmed reads (batch 3) land here later;
    // BWA_ALIGN still consumes EXTRACT_REGION's untrimmed output directly
    // until batch 4 rewires it to the trimmed reads.
    FASTQC('raw', EXTRACT_REGION.out.r1, EXTRACT_REGION.out.r2)
    MULTIQC('raw', FASTQC.out.zip)

    // Batch 2: fastp adapter/quality trimming. Output isn't consumed by
    // BWA_ALIGN yet on purpose — that rewire is batch 4, once batch 3's
    // post-trim FASTQC/MULTIQC pass has confirmed what trimming actually
    // changed. For now TRIM_READS runs alongside the existing pipeline
    // without altering it.
    TRIM_READS(EXTRACT_REGION.out.r1, EXTRACT_REGION.out.r2)

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

    CROSS_CHECK_VCFS(
        GATK_CALL.out.vcf,
        GATK_CALL.out.tbi,
        DEEPVARIANT_CALL.out.vcf,
        DEEPVARIANT_CALL.out.tbi
    )

    FETCH_TRUTH_SET(regions_bed, params.truth_vcf_url, params.truth_bed_url)

    // Truth-set outputs feed two HAPPY_BENCHMARK calls (one per caller), so
    // they need to be value channels for the same reason as
    // regions_bed/reference_fasta/reference_fai above.
    truth_vcf          = FETCH_TRUTH_SET.out.vcf.first()
    truth_tbi          = FETCH_TRUTH_SET.out.tbi.first()
    truth_confident_bed = FETCH_TRUTH_SET.out.confident_bed.first()

    gatk_tuple = GATK_CALL.out.vcf.combine(GATK_CALL.out.tbi).map { vcf, tbi -> tuple('gatk', vcf, tbi) }
    deepvariant_tuple = DEEPVARIANT_CALL.out.vcf.combine(DEEPVARIANT_CALL.out.tbi).map { vcf, tbi -> tuple('deepvariant', vcf, tbi) }

    HAPPY_BENCHMARK(
        gatk_tuple.mix(deepvariant_tuple),
        reference_fasta,
        reference_fai,
        regions_bed,
        truth_vcf,
        truth_tbi,
        truth_confident_bed
    )

    // Annotates the concordant (GATK ∩ DeepVariant) call set — this
    // pipeline's highest-confidence answer — with transcript consequence
    // (VEP) and population frequency (gnomAD v4.1), feeding
    // CAPN3-DMD-variant-classifier's VariantEvidenceBundle adapter.
    gene_annotation_gff     = Channel.value(file(params.gene_annotation_gff))
    gene_annotation_gff_tbi = Channel.value(file("${params.gene_annotation_gff}.tbi"))

    VEP_ANNOTATE(
        CROSS_CHECK_VCFS.out.concordant_gatk_repr,
        gene_annotation_gff,
        gene_annotation_gff_tbi,
        reference_fasta,
        reference_fai
    )

    GNOMAD_ANNOTATE(
        VEP_ANNOTATE.out.vcf,
        VEP_ANNOTATE.out.tbi,
        regions_bed
    )
}
