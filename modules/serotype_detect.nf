process serotype_detect {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/serotype" }, mode: 'copy'

    input:
        tuple val(meta), path(coverage_files)

    output:
        tuple val(meta), path("${meta.id}_serotype.txt"),        emit: serotype
        tuple val(meta), path("${meta.id}_serotype_detail.tsv"), emit: detail

    script:
    def prefix = meta.id
    """
    serotype_detect.py \\
        --denv1        ${prefix}_DENV1.coverage.txt \\
        --denv2        ${prefix}_DENV2.coverage.txt \\
        --denv3        ${prefix}_DENV3.coverage.txt \\
        --denv4        ${prefix}_DENV4.coverage.txt \\
        --sample-id    ${prefix} \\
        --min-coverage 10 \\
        --output       ${prefix}_serotype.txt \\
        --detail       ${prefix}_serotype_detail.tsv
    """
}
