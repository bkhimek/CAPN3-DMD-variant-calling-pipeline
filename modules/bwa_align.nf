process BWA_ALIGN {
    tag "BWA-MEM2 align to GRCh38 chr15+chrX"
    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'
    publishDir "${params.outdir}/bwa_align", mode: 'copy'

    input:
    path r1
    path r2
    path reference_fasta
    path reference_fai
    path reference_0123
    path reference_amb
    path reference_ann
    path reference_bwt
    path reference_pac

    output:
    path "HG002_aligned.sam", emit: sam

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:HG002\\tSM:HG002\\tPL:ILLUMINA\\tLB:HG002_60x' \\
        ${reference_fasta} ${r1} ${r2} > HG002_aligned.sam
    """
}
