process nextclade_download {
    storeDir "${params.assets_dir}/nextclade/${sero}"

    input:
        val sero

    output:
        tuple val(sero), path("nextclade_dataset"), emit: db

    script:
    def name = sero.toLowerCase()
    """
    nextclade dataset get \\
        --name "community/v-gen-lab/dengue/${name}" \\
        --output-dir nextclade_dataset
    """
}

process nextclade {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/nextclade", mode: 'copy'

    input:
        tuple val(meta), path(consensus), path(dataset)

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
