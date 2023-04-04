process fastqc {
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
      echo "serotype 1, ${x}"
      mkdir -p ${params.output}/dengue1
      mkdir -p ${params.output}/dengue1/assemblies
      mkdir -p ${params.output}/dengue1/variants
      mkdir -p ${params.output}/dengue1/vadr_error_reports
      mkdir -p ${params.output}/dengue1/${x}
      cp ${params.input}/${x}_*.fastq.gz ${params.output}/dengue1/${x}
   
      #Run fastqc on original reads
      fastqc ${params.output}/dengue1/${x}/${x}_1.fastq.gz ${params.output}/dengue1/${x}/${x}_2.fastq.gz

      mv ${params.output}/dengue1/${x}/${x}_1_fastqc.html ${params.output}/dengue1/${x}/${x}_1_original_fastqc.html
      mv ${params.output}/dengue1/${x}/${x}_1_fastqc.zip ${params.output}/dengue1/${x}/${x}_1_original_fastqc.zip
      mv ${params.output}/dengue1/${x}/${x}_2_fastqc.html ${params.output}/dengue1/${x}/${x}_2_original_fastqc.html
      mv ${params.output}/dengue1/${x}/${x}_2_fastqc.zip ${params.output}/dengue1/${x}/${x}_2_original_fastqc.zip

      
   elif [[ ${x} =~ SER2_ ]];then
      echo "serotype 2, ${x}"
      mkdir -p ${params.output}/dengue2
      mkdir -p ${params.output}/dengue2/assemblies
      mkdir -p ${params.output}/dengue2/variants
      mkdir -p ${params.output}/dengue2/vadr_error_reports
      mkdir -p ${params.output}/dengue2/${x}
      cp ${params.input}/${x}_*.fastq.gz ${params.output}/dengue2/${x}
   
      #Run fastqc on original reads
      fastqc ${params.output}/dengue2/${x}/${x}_1.fastq.gz ${params.output}/dengue2/${x}/${x}_2.fastq.gz

      mv ${params.output}/dengue2/${x}/${x}_1_fastqc.html ${params.output}/dengue2/${x}/${x}_1_original_fastqc.html
      mv ${params.output}/dengue2/${x}/${x}_1_fastqc.zip ${params.output}/dengue2/${x}/${x}_1_original_fastqc.zip
      mv ${params.output}/dengue2/${x}/${x}_2_fastqc.html ${params.output}/dengue2/${x}/${x}_2_original_fastqc.html
      mv ${params.output}/dengue2/${x}/${x}_2_fastqc.zip ${params.output}/dengue2/${x}/${x}_2_original_fastqc.zip


      
   elif [[ ${x} =~ SER3_ ]];then
      echo "serotype 3, ${x}"
      mkdir -p ${params.output}/dengue3
      mkdir -p ${params.output}/dengue3/assemblies
      mkdir -p ${params.output}/dengue3/variants
      mkdir -p ${params.output}/dengue3/vadr_error_reports
      mkdir -p ${params.output}/dengue3/${x}
      cp ${params.input}/${x}_*.fastq.gz ${params.output}/dengue3/${x}
   
      #Run fastqc on original reads
      fastqc ${params.output}/dengue3/${x}/${x}_1.fastq.gz ${params.output}/dengue3/${x}/${x}_2.fastq.gz


      mv ${params.output}/dengue3/${x}/${x}_1_fastqc.html ${params.output}/dengue3/${x}/${x}_1_original_fastqc.html
      mv ${params.output}/dengue3/${x}/${x}_1_fastqc.zip ${params.output}/dengue3/${x}/${x}_1_original_fastqc.zip
      mv ${params.output}/dengue3/${x}/${x}_2_fastqc.html ${params.output}/dengue3/${x}/${x}_2_original_fastqc.html
      mv ${params.output}/dengue3/${x}/${x}_2_fastqc.zip ${params.output}/dengue3/${x}/${x}_2_original_fastqc.zip
   

      
   elif [[ ${x} =~ SER4_ ]];then
      echo "serotype 4, ${x}"
            mkdir -p ${params.output}/dengue4
      mkdir -p ${params.output}/dengue4/assemblies
      mkdir -p ${params.output}/dengue4/variants
      mkdir -p ${params.output}/dengue4/vadr_error_reports
      mkdir -p ${params.output}/dengue4/${x}
      cp ${params.input}/${x}_*.fastq.gz ${params.output}/dengue4/${x}
   
      #Run fastqc on original reads
      fastqc ${params.output}/dengue4/${x}/${x}_1.fastq.gz ${params.output}/dengue4/${x}/${x}_2.fastq.gz


      mv ${params.output}/dengue4/${x}/${x}_1_fastqc.html ${params.output}/dengue4/${x}/${x}_1_original_fastqc.html
      mv ${params.output}/dengue4/${x}/${x}_1_fastqc.zip ${params.output}/dengue4/${x}/${x}_1_original_fastqc.zip
      mv ${params.output}/dengue4/${x}/${x}_2_fastqc.html ${params.output}/dengue4/${x}/${x}_2_original_fastqc.html
      mv ${params.output}/dengue4/${x}/${x}_2_fastqc.zip ${params.output}/dengue4/${x}/${x}_2_original_fastqc.zip
   

   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
