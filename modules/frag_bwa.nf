process frag_bwa {
    input:
        //val mypath
        val samplename
    output:
        //stdout
        //val mypath
        val samplename
        //path "pyoutputs.txt", emit: pyoutputs

    """
    if [[ ${samplename} =~ SER1_ ]];then
       bwa mem ${params.reference}/NC_001477.1_DENV1.fasta ${params.output}/dengue1/${samplename}/${samplename}_1.fq.gz ${params.output}/dengue1/${samplename}/${samplename}_2.fq.gz > ${params.output}/dengue1/${samplename}/aln-se.sam
       
    elif [[ ${samplename} =~ SER2_ ]];then
       bwa mem ${params.reference}/NC_001474.2_DENV2.fasta ${params.output}/dengue2/${samplename}/${samplename}_1.fq.gz ${params.output}/dengue2/${samplename}/${samplename}_2.fq.gz > ${params.output}/dengue2/${samplename}/aln-se.sam

    elif [[ ${samplename} =~ SER3_ ]];then
       bwa mem ${params.reference}/NC_001475.2_DENV3.fasta ${params.output}/dengue3/${samplename}/${samplename}_1.fq.gz ${params.output}/dengue3/${samplename}/${samplename}_2.fq.gz > ${params.output}/dengue3/${samplename}/aln-se.sam

    elif [[ ${samplename} =~ SER4_ ]];then
       bwa mem ${params.reference}/NC_002640.1_DENV4.fasta ${params.output}/dengue4/${samplename}/${samplename}_1.fq.gz ${params.output}/dengue4/${samplename}/${samplename}_2.fq.gz > ${params.output}/dengue4/${samplename}/aln-se.sam

    else
      echo "check frag.nf step!"
    fi
    
    """
}
