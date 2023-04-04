process primer_trim_samtools {
    input:
        val samplename
 
    output:
        //stdout
        val samplename
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """
    if [[ ${samplename} =~ SER1_ ]];then
       samtools sort ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim.bam -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools index ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools coverage ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim.sorted.bam -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.coverage.txt

       samtools mpileup -A -d 8000 --reference ${params.reference}/NC_001477.1_DENV1.fasta -B -Q 0 -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.mpileup ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim.sorted.bam

    
    elif [[ ${samplename} =~ SER2_ ]];then
       samtools sort ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim.bam -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools index ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools coverage ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim.sorted.bam -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.coverage.txt

       samtools mpileup -A -d 8000 --reference ${params.reference}/NC_001474.2_DENV2.fasta -B -Q 0 -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.mpileup ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim.sorted.bam
    
    elif [[ ${samplename} =~ SER3_ ]];then
       samtools sort ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim.bam -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools index ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools coverage ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim.sorted.bam -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.coverage.txt

       samtools mpileup -A -d 8000 --reference ${params.reference}/NC_001475.2_DENV3.fasta -B -Q 0 -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.mpileup ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim.sorted.bam
    
    elif [[ ${samplename} =~ SER4_ ]];then
       samtools sort ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim.bam -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools index ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim.sorted.bam
       samtools coverage ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim.sorted.bam -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.coverage.txt

       samtools mpileup -A -d 8000 --reference ${params.reference}/NC_002640.1_DENV4.fasta -B -Q 0 -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.mpileup ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim.sorted.bam
    
    else
      echo "check primer.nf step!"
    fi
    
    """
}