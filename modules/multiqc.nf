process multiqc {
    tag "multiqc"
    publishDir "${params.output}/multiqc", mode: 'copy'

    input:
        path(qc_files)
    output:
        path("multiqc_report.html"), emit: report
        path("multiqc_data/"),       emit: data

    script:
    """
    multiqc .
    """
}
