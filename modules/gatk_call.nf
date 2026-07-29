process GATK_CALL {
    tag "GATK HaplotypeCaller"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/gatk_call", mode: 'copy'

    input:
    path bam
    path bai
    path reference_fasta
    path reference_fai
    path reference_dict
    path regions_bed

    output:
    path "HG002_gatk.vcf.gz", emit: vcf
    path "HG002_gatk.vcf.gz.tbi", emit: tbi

    script:
    """
    gatk --java-options "-Xmx6g" HaplotypeCaller \\
        -R ${reference_fasta} \\
        -I ${bam} \\
        -L ${regions_bed} \\
        -O HG002_gatk.vcf.gz
    """
}
