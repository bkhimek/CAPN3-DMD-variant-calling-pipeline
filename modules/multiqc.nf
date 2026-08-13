process MULTIQC {
    tag "MultiQC — ${stage} reads"
    // quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1 — live-verified
    // 2026-08-13 via the quay.io API tag list (1.35 is MultiQC's current
    // release; _1 is the newer of the two 1.35 builds present).
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir { "${params.outdir}/fastq_qc/${stage}/multiqc" }, mode: 'copy'
 
    input:
    val stage
    path fastqc_reports
 
    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
 
    script:
    """
    # MultiQC reads FastQC's zip archives directly (no need to unzip first —
    # it parses the embedded fastqc_data.txt on its own).
    # Note: --filename multiqc_report.html makes MultiQC name its data dir
    # after that filename stem too (multiqc_report_data), not the default
    # multiqc_data — confirmed against actual run output, the output block
    # above must match that or Nextflow fails with "missing output file(s)".
    multiqc . --filename multiqc_report.html
    """
}
