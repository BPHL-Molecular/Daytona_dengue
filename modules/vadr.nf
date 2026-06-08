process vadr {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/vadr", mode: 'copy'

    input:
        tuple val(meta), path(consensus)
    output:
        tuple val(meta), path("vadr_results/"), emit: results
        val meta,                               emit: done

    script:
    def prefix = meta.id
    """
    # Trim terminal ambiguous bases before annotation
    fasta-trim-terminal-ambigs.pl \\
        --minlen 50 \\
        --maxlen 30000 \\
        ${consensus} \\
        > ${prefix}.trimmed.fasta

    # Run VADR annotation
    v-annotate.pl \\
        --split \\
        --cpu ${task.cpus} \\
        --group Dengue \\
        --nomisc \\
        --noprotid \\
        --mkey flavi \\
        --mdir /opt/vadr/vadr-models-flavi/ \\
        --noseqnamemax \\
        ${prefix}.trimmed.fasta \\
        vadr_results
    """
}
