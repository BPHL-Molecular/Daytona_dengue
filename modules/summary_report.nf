process summary_report {
    tag "summary"
    publishDir { "${params.output}" }, mode: 'copy', pattern: 'summary_report.txt'

    input:
        val  barrier
        path qc_files
        path coverage_files
        path consensus_files
        path nextclade_files
        path vadr_dirs
        path kraken2_reports
        path serotype_files
        path screen_cov_files
        path trimstat_files
        path phix_log_files
    output:
        path "summary_report.txt", emit: report
        path "*_mqc.tsv",          emit: mqc_tables, optional: true

    script:
    """
    summary_report.py \\
        --qc-dir         . \\
        --coverage-dir   . \\
        --consensus-dir  . \\
        --nextclade-dir  . \\
        --vadr-dir       . \\
        --kraken2-dir    . \\
        --serotype-dir   . \\
        --screen-cov-dir . \\
        --trimstat-dir   . \\
        --phix-log-dir   . \\
        --output         summary_report.txt
    """
}
