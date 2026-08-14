process TRIM_READS {
    tag "fastp — adapter/quality trimming"
    // quay.io/biocontainers/fastp:1.3.6--h43da1c4_0 — live-verified 2026-08-14
    // directly against the quay.io tag API (matches fastp's actual current
    // release, 1.3.6). Note for the record: an earlier indirect check (a
    // third-party registry mirror + nf-core's own module, both capped at
    // 1.1.0) wrongly suggested no newer build existed and this exact tag
    // was a guess — both of those secondary sources turned out to be
    // stale/lagging, not the tag. The live quay.io API is the actual
    // ground truth here, which is why this project checks it directly
    // rather than trusting indirect sources when they're reachable.
    container 'quay.io/biocontainers/fastp:1.3.6--h43da1c4_0'
    publishDir "${params.outdir}/trim_reads", mode: 'copy'
 
    input:
    path r1
    path r2
 
    output:
    path "HG002_region_trimmed_R1.fastq.gz", emit: r1
    path "HG002_region_trimmed_R2.fastq.gz", emit: r2
    path "fastp_report.html", emit: html
    path "fastp_report.json", emit: json
 
    script:
    """
    # --detect_adapter_for_pe adds fastp's k-mer-based adapter autodetection
    # on top of its default paired-end overlap analysis — belt-and-braces
    # for adapter calling since these reads' true insert size isn't known
    # ahead of time. Everything else left at fastp's defaults (Q15 quality
    # cutoff, standard length filtering) rather than hand-tuned — the point
    # of this batch is to see what fastp's defaults actually do to already
    # high-quality GIAB data, not to force a particular outcome.
    fastp \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        -i ${r1} -I ${r2} \\
        -o HG002_region_trimmed_R1.fastq.gz -O HG002_region_trimmed_R2.fastq.gz \\
        --json fastp_report.json --html fastp_report.html
    """
}
 
