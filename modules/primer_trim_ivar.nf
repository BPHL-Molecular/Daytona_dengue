process primer_trim_ivar {
    input:
        val samplename
 
    output:
        //stdout
        val samplename
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """
    if [[ ${samplename} =~ SER1_ ]];then
       #Trim primers with iVar
       ivar trim -i ${params.output}/dengue1/${samplename}/alignment/${samplename}.sorted.bam  -b ${params.primer}/DENV1.primer.bed -p ${params.output}/dengue1/${samplename}/alignment/${samplename}.primertrim -e

    elif [[ ${samplename} =~ SER2_ ]];then
       #Trim primers with iVar
       ivar trim -i ${params.output}/dengue2/${samplename}/alignment/${samplename}.sorted.bam  -b ${params.primer}/DENV2.primer.bed -p ${params.output}/dengue2/${samplename}/alignment/${samplename}.primertrim -e

    elif [[ ${samplename} =~ SER3_ ]];then
       #Trim primers with iVar
       ivar trim -i ${params.output}/dengue3/${samplename}/alignment/${samplename}.sorted.bam  -b ${params.primer}/DENV3.primer.bed -p ${params.output}/dengue3/${samplename}/alignment/${samplename}.primertrim -e

    elif [[ ${samplename} =~ SER4_ ]];then
       #Trim primers with iVar
       ivar trim -i ${params.output}/dengue4/${samplename}/alignment/${samplename}.sorted.bam  -b ${params.primer}/DENV4.primer.bed -p ${params.output}/dengue4/${samplename}/alignment/${samplename}.primertrim -e

    else
      echo "check primer.nf step!"
    fi
    
    """
}