process EXTRACT_REGION {
    tag "HG002 CAPN3+DMD region extraction"
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    publishDir "${params.outdir}/extract_region", mode: 'copy'

    input:
    path regions_bed
    val bam_url

    output:
    path "HG002_region_R1.fastq.gz", emit: r1
    path "HG002_region_R2.fastq.gz", emit: r2
    path "HG002_region_subset.bam", emit: region_bam

    script:
    """
    # Region strings must be chr-prefixed to match the source BAM's contigs
    # (see docs/regions.md) — samtools does remote random access via the
    # BAM's HTTP-range support and its .bai index, no full download needed.
    # Skip the '#chrom...' header line — samtools silently tolerates it as
    # an invalid-region warning and continues, but that's undocumented
    # leniency, not something to depend on.
    REGIONS=\$(awk '!/^#/ {print \$1":"(\$2+1)"-"\$3}' ${regions_bed} | tr '\\n' ' ')

    samtools view -b -o HG002_region_subset.bam ${bam_url} \$REGIONS
    samtools sort -n -@ ${task.cpus} -o HG002_region_namesorted.bam HG002_region_subset.bam
    samtools fastq -@ ${task.cpus} \\
        -1 HG002_region_R1.fastq.gz \\
        -2 HG002_region_R2.fastq.gz \\
        -0 /dev/null -s /dev/null -n \\
        HG002_region_namesorted.bam
    """
}
