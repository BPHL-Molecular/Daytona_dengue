process multiqc {
   input:
      //tuple val(y), val(x)
      val(x)

   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      //val "${params.output}/dengue1/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      val "${x}"
      
   """  
   if [[ ${x} =~ SER1_ ]];then
      #Run multiqc
      multiqc ${params.output}/dengue1/${x}/${x}_*_fastqc.zip -o ${params.output}/dengue1/${x}

      #Map reads to reference
      mkdir ${params.output}/dengue1/${x}/alignment
      mv ${params.output}/dengue1/${x}/multiqc_data ${params.output}/dengue1/${x}/${x}_multiqc_data 
      
   elif [[ ${x} =~ SER2_ ]];then
      #Run multiqc
      multiqc ${params.output}/dengue2/${x}/${x}_*_fastqc.zip -o ${params.output}/dengue2/${x}

      #Map reads to reference
      mkdir ${params.output}/dengue2/${x}/alignment
      mv ${params.output}/dengue2/${x}/multiqc_data ${params.output}/dengue2/${x}/${x}_multiqc_data
      
   elif [[ ${x} =~ SER3_ ]];then
      #Run multiqc
      multiqc ${params.output}/dengue3/${x}/${x}_*_fastqc.zip -o ${params.output}/dengue3/${x}

      #Map reads to reference
      mkdir ${params.output}/dengue3/${x}/alignment
      mv ${params.output}/dengue3/${x}/multiqc_data ${params.output}/dengue3/${x}/${x}_multiqc_data
      
   elif [[ ${x} =~ SER4_ ]];then
      #Run multiqc
      multiqc ${params.output}/dengue4/${x}/${x}_*_fastqc.zip -o ${params.output}/dengue4/${x}

      #Map reads to reference
      mkdir ${params.output}/dengue4/${x}/alignment
      mv ${params.output}/dengue4/${x}/multiqc_data ${params.output}/dengue4/${x}/${x}_multiqc_data
   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
