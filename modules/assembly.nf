process assembly {
    input:
        val samplename
    output:
        //stdout
        val samplename
        //path "pyoutputs.txt", emit: pyoutputs
        
    
    """ 
    if [[ ${samplename} =~ SER1_ ]];then
        #Call variants
        mkdir ${params.output}/dengue1/${samplename}/variants
        cat ${params.output}/dengue1/${samplename}/alignment/${samplename}.mpileup | tee | ivar variants -r ${params.reference}/NC_001477.1_DENV1.fasta -m 10 -p ${params.output}/dengue1/${samplename}/variants/${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/Dengue1_GCF_000862125.1_ViralProj15306_genomic.gff.gz

        #Generate consensus assembly
        mkdir ${params.output}/dengue1/${samplename}/assembly
        cat ${params.output}/dengue1/${samplename}/alignment/${samplename}.mpileup | tee | ivar consensus -t 0 -m 10 -n N -p ${params.output}/dengue1/${samplename}/assembly/${samplename}.consensus

    elif [[ ${samplename} =~ SER2_ ]];then
        #Call variants
        mkdir ${params.output}/dengue2/${samplename}/variants
        cat ${params.output}/dengue2/${samplename}/alignment/${samplename}.mpileup | tee | ivar variants -r ${params.reference}/NC_001474.2_DENV2.fasta -m 10 -p ${params.output}/dengue2/${samplename}/variants/${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/Dengue2_GCF_000871845.1_ViralProj20183_genomic.gff.gz
        #Generate consensus assembly
        mkdir ${params.output}/dengue2/${samplename}/assembly
        cat ${params.output}/dengue2/${samplename}/alignment/${samplename}.mpileup | tee | ivar consensus -t 0 -m 10 -n N -p ${params.output}/dengue2/${samplename}/assembly/${samplename}.consensus
    
    elif [[ ${samplename} =~ SER3_ ]];then
        #Call variants
        mkdir ${params.output}/dengue3/${samplename}/variants
        cat ${params.output}/dengue3/${samplename}/alignment/${samplename}.mpileup | tee | ivar variants -r ${params.reference}/NC_001475.2_DENV3.fasta -m 10 -p ${params.output}/dengue3/${samplename}/variants/${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/Dengue3_GCF_000866625.1_ViralProj15598_genomic.gff.gz
        #Generate consensus assembly
        mkdir ${params.output}/dengue3/${samplename}/assembly
        cat ${params.output}/dengue3/${samplename}/alignment/${samplename}.mpileup | tee | ivar consensus -t 0 -m 10 -n N -p ${params.output}/dengue3/${samplename}/assembly/${samplename}.consensus
    
    elif [[ ${samplename} =~ SER4_ ]];then
        #Call variants
        mkdir ${params.output}/dengue4/${samplename}/variants
        cat ${params.output}/dengue4/${samplename}/alignment/${samplename}.mpileup | tee | ivar variants -r ${params.reference}/NC_002640.1_DENV4.fasta -m 10 -p ${params.output}/dengue4/${samplename}/variants/${samplename}.variants -q 20 -t 0.25 -g ${params.reference}/Dengue4_GCF_000865065.1_ViralProj15599_genomic.gff.gz
        #Generate consensus assembly
        mkdir ${params.output}/dengue4/${samplename}/assembly
        cat ${params.output}/dengue4/${samplename}/alignment/${samplename}.mpileup | tee | ivar consensus -t 0 -m 10 -n N -p ${params.output}/dengue4/${samplename}/assembly/${samplename}.consensus
    
    else
      echo "check assembly.nf step!"
    fi
    
    
    
    """
}