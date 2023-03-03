process primer2 {
    input:
        val mypath
 
    output:
        //stdout
        val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #Index final sorted bam
    singularity exec docker://staphb/samtools:1.12 samtools index ${mypath}/alignment/\${samplename}.sorted.bam

    
    #Trim primers with iVar
    singularity exec docker://staphb/ivar ivar trim -i ${mypath}/alignment/\${samplename}.sorted.bam  -b ${params.primer}/DENV2.primer.bed -p ${mypath}/alignment/\${samplename}.primertrim -e
    
    singularity exec docker://staphb/samtools:1.12 samtools sort ${mypath}/alignment/\${samplename}.primertrim.bam -o ${mypath}/alignment/\${samplename}.primertrim.sorted.bam

    singularity exec docker://staphb/samtools:1.12 samtools index ${mypath}/alignment/\${samplename}.primertrim.sorted.bam

    singularity exec docker://staphb/samtools:1.12 samtools coverage ${mypath}/alignment/\${samplename}.primertrim.sorted.bam -o ${mypath}/alignment/\${samplename}.coverage.txt
    """
}