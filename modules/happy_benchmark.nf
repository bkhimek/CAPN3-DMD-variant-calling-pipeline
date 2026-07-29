process HAPPY_BENCHMARK {
    tag "hap.py benchmark: ${caller}"
    container 'jmcdani20/hap.py:v0.3.12'
    publishDir "${params.outdir}/happy_benchmark", mode: 'copy', saveAs: { fn -> "${caller}/${fn}" }

    input:
    tuple val(caller), path(query_vcf), path(query_tbi)
    path reference_fasta
    path reference_fai
    path regions_bed
    path truth_vcf
    path truth_tbi
    path truth_confident_bed

    output:
    path "${caller}.summary.csv", emit: summary
    path "${caller}.extended.csv", emit: extended
    path "${caller}.vcf.gz", emit: annotated_vcf

    script:
    """
    # --pass-only verified (interactive pretest, batch 7) to correctly keep
    # GATK's FILTER='.' (unfiltered/final) records rather than wrongly
    # excluding them — no FILTER normalization needed.
    # -f = confident/callable regions (where truth is trustworthy); -T = our
    # padded CAPN3/DMD target regions (dense — a couple of large contiguous
    # regions — hence -T rather than -R, per hap.py's own guidance).
    /opt/hap.py/bin/hap.py \\
        ${truth_vcf} ${query_vcf} \\
        -r ${reference_fasta} \\
        -f ${truth_confident_bed} \\
        -T ${regions_bed} \\
        -o ${caller} \\
        --pass-only
    """
}
