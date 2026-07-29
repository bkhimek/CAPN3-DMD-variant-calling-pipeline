process SORT_MARKDUP {
    tag "coordinate sort + mark duplicates"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/sort_markdup", mode: 'copy'

    input:
    path sam

    output:
    path "HG002_markdup.bam", emit: bam
    path "HG002_markdup.bam.bai", emit: bai
    path "HG002_markdup_metrics.txt", emit: metrics

    script:
    """
    # broadinstitute/gatk bundles samtools too, so sort/index and MarkDuplicates
    # can share one container — no separate Picard-only image needed.
    samtools sort -@ ${task.cpus} -o HG002_sorted.bam ${sam}

    gatk --java-options "-Xmx6g" MarkDuplicates \\
        -I HG002_sorted.bam \\
        -O HG002_markdup.bam \\
        -M HG002_markdup_metrics.txt

    samtools index HG002_markdup.bam
    """
}
