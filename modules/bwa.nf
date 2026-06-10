process bwa {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/bwa", mode: 'copy', pattern: "*.sam"

    input:
        tuple val(meta), path(reads)
        path denv1_files
        path denv2_files
        path denv3_files
        path denv4_files

    output:
        tuple val(meta), path("*_DENV?.sam"), emit: sams

    script:
    def prefix   = meta.id
    def denv1_fa = denv1_files instanceof List ? denv1_files.find { f -> f.name.endsWith('.fasta') } : denv1_files
    def denv2_fa = denv2_files instanceof List ? denv2_files.find { f -> f.name.endsWith('.fasta') } : denv2_files
    def denv3_fa = denv3_files instanceof List ? denv3_files.find { f -> f.name.endsWith('.fasta') } : denv3_files
    def denv4_fa = denv4_files instanceof List ? denv4_files.find { f -> f.name.endsWith('.fasta') } : denv4_files
    """
    bwa mem -t ${task.cpus} ${denv1_fa} ${reads[0]} ${reads[1]} > ${prefix}_DENV1.sam
    bwa mem -t ${task.cpus} ${denv2_fa} ${reads[0]} ${reads[1]} > ${prefix}_DENV2.sam
    bwa mem -t ${task.cpus} ${denv3_fa} ${reads[0]} ${reads[1]} > ${prefix}_DENV3.sam
    bwa mem -t ${task.cpus} ${denv4_fa} ${reads[0]} ${reads[1]} > ${prefix}_DENV4.sam
    """
}
