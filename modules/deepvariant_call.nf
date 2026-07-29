process DEEPVARIANT_CALL {
    tag "DeepVariant WGS caller"
    container 'google/deepvariant:1.10.0'
    publishDir "${params.outdir}/deepvariant_call", mode: 'copy'

    input:
    path bam
    path bai
    path reference_fasta
    path reference_fai
    path regions_bed

    output:
    path "HG002_deepvariant.vcf.gz", emit: vcf
    path "HG002_deepvariant.vcf.gz.tbi", emit: tbi

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${reference_fasta} \\
        --reads=${bam} \\
        --regions=${regions_bed} \\
        --output_vcf=HG002_deepvariant.vcf.gz \\
        --num_shards=${task.cpus}
    """
}
