process FASTQC {
    tag "FastQC — ${stage} reads"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir { "${params.outdir}/fastq_qc/${stage}/fastqc" }, mode: 'copy'

    input:
    val stage
    path r1
    path r2

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    # Runs against R1 and R2 in one invocation — FastQC reports per-file, so
    # this always produces two zip/html pairs (one per read direction), not
    # one merged report. `stage` (raw vs trimmed) only controls where the
    # output lands via publishDir; it isn't passed to fastqc itself.
    fastqc --threads ${task.cpus} ${r1} ${r2}
    """
}
