process frag2 {
    input:
        val mypath
 
    output:
        //stdout
        val mypath
        //path "pyoutputs.txt", emit: pyoutputs

    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #singularity exec docker://staphb/bwa:0.7.17 bwa index ${params.reference}/NC_001474.2_DENV2.fasta
    singularity exec docker://staphb/bwa:0.7.17 bwa mem ${params.reference}/NC_001474.2_DENV2.fasta ${mypath}/\${samplename}_1.fq.gz ${mypath}/\${samplename}_2.fq.gz | singularity exec docker://staphb/samtools:1.12 samtools view - -F 4 -u -h | singularity exec docker://staphb/samtools:1.12 samtools sort -n > ${mypath}/alignment/\${samplename}.namesorted.bam

    singularity exec docker://staphb/samtools:1.12 samtools fixmate -m ${mypath}/alignment/\${samplename}.namesorted.bam ${mypath}/alignment/\${samplename}.fixmate.bam

        #Create positional sorted bam from fixmate.bam
    singularity exec docker://staphb/samtools:1.12 samtools sort -o ${mypath}/alignment/\${samplename}.positionsort.bam ${mypath}/alignment/\${samplename}.fixmate.bam

        #Mark duplicate reads
    singularity exec docker://staphb/samtools:1.12 samtools markdup ${mypath}/alignment/\${samplename}.positionsort.bam ${mypath}/alignment/\${samplename}.markdup.bam

        #Remove duplicate reads
    singularity exec docker://staphb/samtools:1.12 samtools markdup -r ${mypath}/alignment/\${samplename}.positionsort.bam ${mypath}/alignment/\${samplename}.dedup.bam

        #Sort dedup.bam and rename to .sorted.bam
    singularity exec docker://staphb/samtools:1.12 samtools sort -o ${mypath}/alignment/\${samplename}.sorted.bam ${mypath}/alignment/\${samplename}.dedup.bam

       #Index final sorted bam
    #singularity exec /apps/staphb-toolkit/containers/samtools_1.12.sif samtools index ${mypath}/alignment/\${samplename}.sorted.bam
    """
}