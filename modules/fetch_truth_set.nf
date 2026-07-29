process FETCH_TRUTH_SET {
    tag "GIAB HG002 truth set — CAPN3+DMD padded regions"
    container 'quay.io/biocontainers/bcftools:1.24--h487d631_1'
    publishDir "${params.outdir}/happy_benchmark", mode: 'copy'

    input:
    path regions_bed
    val truth_vcf_url
    val truth_bed_url

    output:
    path "truth_region.vcf.gz", emit: vcf
    path "truth_region.vcf.gz.tbi", emit: tbi
    path "truth_confident.bed", emit: confident_bed

    script:
    """
    # Truth VCF has a .tbi and supports HTTP byte ranges (same trick as the
    # HG002 BAM in EXTRACT_REGION) — pull just the padded regions instead of
    # the full 156MB genome-wide file.
    REGIONS=\$(awk '!/^#/ {print \$1":"(\$2+1)"-"\$3}' ${regions_bed} | tr '\\n' ' ')
    bcftools view -Oz -o truth_region.vcf.gz ${truth_vcf_url} \$REGIONS
    bcftools index -t truth_region.vcf.gz

    # Confident-regions BED is plain (uncompressed) text, no random-access
    # option, but at ~11MB it's cheap to pull in full.
    wget -q -O truth_confident.bed ${truth_bed_url}
    """
}
