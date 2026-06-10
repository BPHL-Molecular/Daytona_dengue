process vadr_download {
    storeDir "${params.assets_dir}/vadr"

    output:
        path "vadr-models-flavi-1.2-1", emit: models

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flaviviridae/CURRENT/vadr-models-flavi-1.2-1.tar.gz
    tar -xf vadr-models-flavi-1.2-1.tar.gz
    rm vadr-models-flavi-1.2-1.tar.gz
    """
}

process vadr {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/vadr", mode: 'copy'

    input:
        tuple val(meta), path(consensus)
        path(vadr_models)
    output:
        tuple val(meta), path("${meta.id}_vadr_results/"), emit: results
        val meta,                                          emit: done

    script:
    def prefix = meta.id
    """
    fasta-trim-terminal-ambigs.pl \\
        --minlen 50 \\
        --maxlen 30000 \\
        ${consensus} \\
        > ${prefix}.trimmed.fasta

    v-annotate.pl \\
        --split \\
        --cpu ${task.cpus} \\
        --group Dengue \\
        --nomisc \\
        --noprotid \\
        --mkey flavi \\
        --mdir ${vadr_models} \\
        --noseqnamemax \\
        ${prefix}.trimmed.fasta \\
        ${prefix}_vadr_results
    """
}
