process serotype_detect {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/serotype", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        path denv1_files   // [fasta, .amb, .ann, .bwt, .fai, .pac, .sa] staged flat
        path denv2_files
        path denv3_files
        path denv4_files

    output:
        tuple val(meta), path("${meta.id}_serotype.txt"),        emit: serotype
        tuple val(meta), path("${meta.id}_serotype_detail.tsv"), emit: detail

    script:
    def prefix   = meta.id
    // When multiple files are staged, Nextflow presents them as a list;
    // the FASTA is first (alphabetically, .fasta sorts before index extensions)
    def denv1_fa = denv1_files instanceof List ? denv1_files.find { f -> f.name.endsWith('.fasta') } : denv1_files
    def denv2_fa = denv2_files instanceof List ? denv2_files.find { f -> f.name.endsWith('.fasta') } : denv2_files
    def denv3_fa = denv3_files instanceof List ? denv3_files.find { f -> f.name.endsWith('.fasta') } : denv3_files
    def denv4_fa = denv4_files instanceof List ? denv4_files.find { f -> f.name.endsWith('.fasta') } : denv4_files
    """
    bwa mem -t ${task.cpus} ${denv1_fa} ${reads[0]} ${reads[1]} \\
        | samtools view -F 4 -b \\
        | samtools coverage -o ${prefix}_DENV1.coverage.txt

    bwa mem -t ${task.cpus} ${denv2_fa} ${reads[0]} ${reads[1]} \\
        | samtools view -F 4 -b \\
        | samtools coverage -o ${prefix}_DENV2.coverage.txt

    bwa mem -t ${task.cpus} ${denv3_fa} ${reads[0]} ${reads[1]} \\
        | samtools view -F 4 -b \\
        | samtools coverage -o ${prefix}_DENV3.coverage.txt

    bwa mem -t ${task.cpus} ${denv4_fa} ${reads[0]} ${reads[1]} \\
        | samtools view -F 4 -b \\
        | samtools coverage -o ${prefix}_DENV4.coverage.txt

    serotype_detect.py \\
        --denv1        ${prefix}_DENV1.coverage.txt \\
        --denv2        ${prefix}_DENV2.coverage.txt \\
        --denv3        ${prefix}_DENV3.coverage.txt \\
        --denv4        ${prefix}_DENV4.coverage.txt \\
        --sample-id    ${prefix} \\
        --min-coverage 10 \\
        --output       ${meta.id}_serotype.txt \\
        --detail       ${meta.id}_serotype_detail.tsv
    """
}
