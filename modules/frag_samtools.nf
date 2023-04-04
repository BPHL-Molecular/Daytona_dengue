process frag_samtools {
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

       samtools view -F 4 -u -h -bo ${params.output}/dengue1/${samplename}/aln-se.bam ${params.output}/dengue1/${samplename}/aln-se.sam  
       samtools sort -n -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue1/${samplename}/aln-se.bam
       
       samtools fixmate -m ${params.output}/dengue1/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue1/${samplename}/alignment/${samplename}.fixmate.bam

        #Create positional sorted bam from fixmate.bam
       samtools sort -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue1/${samplename}/alignment/${samplename}.fixmate.bam

        #Mark duplicate reads
       samtools markdup ${params.output}/dengue1/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue1/${samplename}/alignment/${samplename}.markdup.bam

        #Remove duplicate reads
       samtools markdup -r ${params.output}/dengue1/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue1/${samplename}/alignment/${samplename}.dedup.bam

        #Sort dedup.bam and rename to .sorted.bam
       samtools sort -o ${params.output}/dengue1/${samplename}/alignment/${samplename}.sorted.bam ${params.output}/dengue1/${samplename}/alignment/${samplename}.dedup.bam

        #Index final sorted bam
       samtools index ${params.output}/dengue1/${samplename}/alignment/${samplename}.sorted.bam

       
    elif [[ ${samplename} =~ SER2_ ]];then
       samtools view -F 4 -u -h -bo ${params.output}/dengue2/${samplename}/aln-se.bam ${params.output}/dengue2/${samplename}/aln-se.sam  
       samtools sort -n -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue2/${samplename}/aln-se.bam

       
       samtools fixmate -m ${params.output}/dengue2/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue2/${samplename}/alignment/${samplename}.fixmate.bam

        #Create positional sorted bam from fixmate.bam
       samtools sort -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue2/${samplename}/alignment/${samplename}.fixmate.bam

        #Mark duplicate reads
       samtools markdup ${params.output}/dengue2/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue2/${samplename}/alignment/${samplename}.markdup.bam

        #Remove duplicate reads
       samtools markdup -r ${params.output}/dengue2/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue2/${samplename}/alignment/${samplename}.dedup.bam

        #Sort dedup.bam and rename to .sorted.bam
       samtools sort -o ${params.output}/dengue2/${samplename}/alignment/${samplename}.sorted.bam ${params.output}/dengue2/${samplename}/alignment/${samplename}.dedup.bam

        #Index final sorted bam
       samtools index ${params.output}/dengue2/${samplename}/alignment/${samplename}.sorted.bam

    elif [[ ${samplename} =~ SER3_ ]];then
       samtools view -F 4 -u -h -bo ${params.output}/dengue3/${samplename}/aln-se.bam ${params.output}/dengue3/${samplename}/aln-se.sam  
       samtools sort -n -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue3/${samplename}/aln-se.bam
       
       samtools fixmate -m ${params.output}/dengue3/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue3/${samplename}/alignment/${samplename}.fixmate.bam

        #Create positional sorted bam from fixmate.bam
       samtools sort -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue3/${samplename}/alignment/${samplename}.fixmate.bam

        #Mark duplicate reads
       samtools markdup ${params.output}/dengue3/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue3/${samplename}/alignment/${samplename}.markdup.bam

        #Remove duplicate reads
       samtools markdup -r ${params.output}/dengue3/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue3/${samplename}/alignment/${samplename}.dedup.bam

        #Sort dedup.bam and rename to .sorted.bam
       samtools sort -o ${params.output}/dengue3/${samplename}/alignment/${samplename}.sorted.bam ${params.output}/dengue3/${samplename}/alignment/${samplename}.dedup.bam

        #Index final sorted bam
       samtools index ${params.output}/dengue3/${samplename}/alignment/${samplename}.sorted.bam

    elif [[ ${samplename} =~ SER4_ ]];then
       samtools view -F 4 -u -h -bo ${params.output}/dengue4/${samplename}/aln-se.bam ${params.output}/dengue4/${samplename}/aln-se.sam  
       samtools sort -n -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue4/${samplename}/aln-se.bam

       samtools fixmate -m ${params.output}/dengue4/${samplename}/alignment/${samplename}.namesorted.bam ${params.output}/dengue4/${samplename}/alignment/${samplename}.fixmate.bam

        #Create positional sorted bam from fixmate.bam
       samtools sort -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue4/${samplename}/alignment/${samplename}.fixmate.bam

        #Mark duplicate reads
       samtools markdup ${params.output}/dengue4/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue4/${samplename}/alignment/${samplename}.markdup.bam

        #Remove duplicate reads
       samtools markdup -r ${params.output}/dengue4/${samplename}/alignment/${samplename}.positionsort.bam ${params.output}/dengue4/${samplename}/alignment/${samplename}.dedup.bam

        #Sort dedup.bam and rename to .sorted.bam
       samtools sort -o ${params.output}/dengue4/${samplename}/alignment/${samplename}.sorted.bam ${params.output}/dengue4/${samplename}/alignment/${samplename}.dedup.bam

        #Index final sorted bam
       samtools index ${params.output}/dengue4/${samplename}/alignment/${samplename}.sorted.bam

    else
      echo "check frag.nf step!"
    fi
    
    """
}
