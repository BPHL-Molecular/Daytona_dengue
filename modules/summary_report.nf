process summary_report {
    tag "summary"
    publishDir "${params.output}", mode: 'copy'

    input:
        val  barrier
        path qc_files
        path coverage_files
        path consensus_files
        path nextclade_files
        path vadr_dirs
        path kraken2_reports
    output:
        path "summary_report.tsv", emit: report

    script:
    """
    summary_report.py \\
        --qc-dir        . \\
        --coverage-dir  . \\
        --consensus-dir . \\
        --nextclade-dir . \\
        --vadr-dir      . \\
        --kraken2-dir   . \\
        --output        summary_report.tsv
    """
}
