process CROSS_CHECK_VCFS {
    tag "GATK vs DeepVariant concordance (bcftools isec)"
    container 'quay.io/biocontainers/bcftools:1.24--h487d631_1'
    publishDir "${params.outdir}/cross_check_vcfs", mode: 'copy'

    input:
    path gatk_vcf
    path gatk_tbi
    path deepvariant_vcf
    path deepvariant_tbi

    output:
    path "isec/0000.vcf", emit: gatk_only
    path "isec/0001.vcf", emit: deepvariant_only
    path "isec/0002.vcf", emit: concordant_gatk_repr
    path "isec/0003.vcf", emit: concordant_deepvariant_repr
    path "concordance_summary.txt", emit: summary

    script:
    """
    # bcftools isec compares by position+REF+ALT: 0000/0001 are the
    # records private to each caller (discordant), 0002/0003 are the shared
    # (concordant) records, once using each caller's own record formatting.
    bcftools isec -p isec ${gatk_vcf} ${deepvariant_vcf}

    DV_ONLY_TOTAL=\$(grep -vc '^#' isec/0001.vcf)
    # DeepVariant emits a record (RefCall/NoCall) for every candidate site it
    # examines, not just confident variant calls, so a large share of its
    # "private" records are non-PASS sites rather than genuine competing
    # variant calls against GATK. Break that out explicitly instead of
    # reporting the raw isec count as if it were pure caller disagreement.
    DV_ONLY_PASS=\$(grep -v '^#' isec/0001.vcf | awk -F'\\t' '\$7 == "PASS"' | wc -l)

    {
        echo "GATK_CALL total:                    \$(bcftools view -H ${gatk_vcf} | wc -l)"
        echo "DEEPVARIANT_CALL total:              \$(bcftools view -H ${deepvariant_vcf} | wc -l)"
        echo "GATK-only (discordant):              \$(grep -vc '^#' isec/0000.vcf)"
        echo "DeepVariant-only, raw (discordant):  \${DV_ONLY_TOTAL}"
        echo "DeepVariant-only, PASS-filtered only: \${DV_ONLY_PASS}"
        echo "  (remainder of DeepVariant-only records are RefCall/NoCall —"
        echo "  DeepVariant reports every candidate site, not just variant"
        echo "  calls, so most of the raw DeepVariant-only count isn't"
        echo "  genuine caller disagreement)"
        echo "Concordant (shared):                 \$(grep -vc '^#' isec/0002.vcf)"
    } > concordance_summary.txt
    """
}
