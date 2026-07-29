process VEP_ANNOTATE {
    tag "VEP transcript consequence annotation"
    container 'ensemblorg/ensembl-vep:release_116.0'
    // Must run as the invoking host uid/gid: the VEP image's default user
    // (uid 999) has no write access to a bind-mounted host directory owned
    // by the Nextflow work-dir's uid, and fails opening its own warnings
    // file otherwise (confirmed live during this module's pretest).
    containerOptions "--user \$(id -u):\$(id -g)"
    publishDir "${params.outdir}/annotate_calls", mode: 'copy'

    input:
    path concordant_vcf
    path gene_annotation_gff
    path gene_annotation_gff_tbi
    path reference_fasta
    path reference_fai

    output:
    path "vep_annotated.vcf.gz", emit: vcf
    path "vep_annotated.vcf.gz.tbi", emit: tbi

    script:
    """
    bgzip -c ${concordant_vcf} > concordant.vcf.gz
    tabix -p vcf concordant.vcf.gz

    vep --input_file concordant.vcf.gz --format vcf \\
        --gff ${gene_annotation_gff} --fasta ${reference_fasta} \\
        --vcf --force_overwrite --no_stats --hgvs --symbol \\
        --output_file vep_annotated_unzipped.vcf

    bgzip -c vep_annotated_unzipped.vcf > vep_annotated.vcf.gz
    tabix -p vcf vep_annotated.vcf.gz
    """
}
