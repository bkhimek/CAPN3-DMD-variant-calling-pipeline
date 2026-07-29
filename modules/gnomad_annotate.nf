process GNOMAD_ANNOTATE {
    tag "gnomAD v4.1 population frequency annotation"
    container 'quay.io/biocontainers/bcftools:1.24--h487d631_1'
    publishDir "${params.outdir}/annotate_calls", mode: 'copy'

    input:
    path vep_vcf
    path vep_tbi
    path regions_bed

    output:
    path "annotated_calls.vcf.gz", emit: vcf
    path "annotated_calls.vcf.gz.tbi", emit: tbi

    script:
    """
    # One remote region query per padded gene region (2 total) rather than
    # per-variant lookups — gnomAD's public sites VCFs support HTTP range +
    # tabix random access (confirmed live), same trick already proven for
    # the HG002 BAM and GIAB truth VCF elsewhere in this pipeline.
    CAPN3_REGION="chr15:\$(awk -F'\\t' '\$1=="chr15"{print \$2+1"-"\$3}' ${regions_bed})"
    DMD_REGION="chrX:\$(awk -F'\\t' '\$1=="chrX"{print \$2+1"-"\$3}' ${regions_bed})"

    bcftools view -Oz -o gnomad_capn3.vcf.gz "${params.gnomad_chr15_url}" "\${CAPN3_REGION}"
    bcftools view -Oz -o gnomad_dmd.vcf.gz "${params.gnomad_chrx_url}" "\${DMD_REGION}"
    tabix -p vcf gnomad_capn3.vcf.gz
    tabix -p vcf gnomad_dmd.vcf.gz
    bcftools concat -a gnomad_capn3.vcf.gz gnomad_dmd.vcf.gz -Oz -o gnomad_region.vcf.gz
    tabix -p vcf gnomad_region.vcf.gz

    # Renamed with a gnomAD_ prefix (source file's own AC/AN/AF/etc. would
    # otherwise collide with GATK's identically-named INFO fields already on
    # these records). AF_grpmax = gnomAD v4's "max AF across ancestry
    # groups" field (the modern name for the old popmax AF) — confirmed via
    # the source VCF's own header, not assumed from memory.
    bcftools annotate -a gnomad_region.vcf.gz \\
        -c 'INFO/gnomAD_AC:=INFO/AC,INFO/gnomAD_AN:=INFO/AN,INFO/gnomAD_AF:=INFO/AF,INFO/gnomAD_AF_grpmax:=INFO/AF_grpmax,INFO/gnomAD_nhomalt:=INFO/nhomalt' \\
        -Oz -o annotated_calls.vcf.gz \\
        ${vep_vcf}
    tabix -p vcf annotated_calls.vcf.gz
    """
}
