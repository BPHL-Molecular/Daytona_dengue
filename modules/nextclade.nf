process nextclade_download {
    storeDir "${params.output}/db/nextclade"

    output:
        path "nextclade_dataset", emit: db

    script:
    """
    nextclade dataset get \\
        --name 'nextstrain/dengue/all' \\
        --output-dir nextclade_dataset
    """
}

process nextclade {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/nextclade", mode: 'copy'

    input:
        tuple val(meta), path(consensus)
        path(dataset)
    output:
        tuple val(meta), path("${meta.id}_nextclade.tsv"), emit: tsv
        val meta,                                          emit: done

    script:
    def prefix = meta.id
    """
    nextclade run \\
        --input-dataset ${dataset} \\
        --output-tsv ${prefix}_nextclade.tsv \\
        ${consensus}
    """
}
