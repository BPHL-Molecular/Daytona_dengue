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
    publishDir { "${params.output}/${meta.id}/vadr" },    mode: 'copy', pattern: "*_vadr_results"
    publishDir { "${params.output}/assemblies_qc_pass" }, mode: 'copy', pattern: "{pass,review}/*.consensus.fasta"

    input:
        tuple val(meta), path(consensus)
        path(vadr_models)
    output:
        tuple val(meta), path("${meta.id}_vadr_results/"),     emit: results
        path "pass/${meta.id}.consensus.fasta",   optional: true
        path "review/${meta.id}.consensus.fasta", optional: true
        val meta,                                              emit: done

    script:
    def prefix = meta.id
    def review_min_cov = 80
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

    pass_fa=\$(find ${prefix}_vadr_results -name '*.vadr.pass.fa' | head -n 1)
    fail_fa=\$(find ${prefix}_vadr_results -name '*.vadr.fail.fa' | head -n 1)
    if [ -s "\$pass_fa" ]; then
        mkdir -p pass && cp "\$pass_fa" pass/${prefix}.consensus.fasta
    elif [ -s "\$fail_fa" ] && ${meta.pct_assembled >= review_min_cov}; then
        mkdir -p review && cp "\$fail_fa" review/${prefix}.consensus.fasta
    fi
    """
}
