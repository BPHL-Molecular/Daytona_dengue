process fastqc_clean {
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
      #Run fastqc on clean forward and reverse reads
      fastqc ${params.output}/dengue1/${x}/${x}_1.fq.gz ${params.output}/dengue1/${x}/${x}_2.fq.gz

      #Rename fastqc output files
      mv ${params.output}/dengue1/${x}/${x}_1_fastqc.html ${params.output}/dengue1/${x}/${x}_1_clean_fastqc.html
      mv ${params.output}/dengue1/${x}/${x}_1_fastqc.zip ${params.output}/dengue1/${x}/${x}_1_clean_fastqc.zip
      mv ${params.output}/dengue1/${x}/${x}_2_fastqc.html ${params.output}/dengue1/${x}/${x}_2_clean_fastqc.html
      mv ${params.output}/dengue1/${x}/${x}_2_fastqc.zip ${params.output}/dengue1/${x}/${x}_2_clean_fastqc.zip


      
   elif [[ ${x} =~ SER2_ ]];then
      #Run fastqc on clean forward and reverse reads
      fastqc ${params.output}/dengue2/${x}/${x}_1.fq.gz ${params.output}/dengue2/${x}/${x}_2.fq.gz

      #Rename fastqc output files
      mv ${params.output}/dengue2/${x}/${x}_1_fastqc.html ${params.output}/dengue2/${x}/${x}_1_clean_fastqc.html
      mv ${params.output}/dengue2/${x}/${x}_1_fastqc.zip ${params.output}/dengue2/${x}/${x}_1_clean_fastqc.zip
      mv ${params.output}/dengue2/${x}/${x}_2_fastqc.html ${params.output}/dengue2/${x}/${x}_2_clean_fastqc.html
      mv ${params.output}/dengue2/${x}/${x}_2_fastqc.zip ${params.output}/dengue2/${x}/${x}_2_clean_fastqc.zip



      
   elif [[ ${x} =~ SER3_ ]];then
      #Run fastqc on clean forward and reverse reads
      fastqc ${params.output}/dengue3/${x}/${x}_1.fq.gz ${params.output}/dengue3/${x}/${x}_2.fq.gz

      #Rename fastqc output files
      mv ${params.output}/dengue3/${x}/${x}_1_fastqc.html ${params.output}/dengue3/${x}/${x}_1_clean_fastqc.html
      mv ${params.output}/dengue3/${x}/${x}_1_fastqc.zip ${params.output}/dengue3/${x}/${x}_1_clean_fastqc.zip
      mv ${params.output}/dengue3/${x}/${x}_2_fastqc.html ${params.output}/dengue3/${x}/${x}_2_clean_fastqc.html
      mv ${params.output}/dengue3/${x}/${x}_2_fastqc.zip ${params.output}/dengue3/${x}/${x}_2_clean_fastqc.zip

   

      
   elif [[ ${x} =~ SER4_ ]];then
      #Run fastqc on clean forward and reverse reads
      fastqc ${params.output}/dengue4/${x}/${x}_1.fq.gz ${params.output}/dengue4/${x}/${x}_2.fq.gz

      #Rename fastqc output files
      mv ${params.output}/dengue4/${x}/${x}_1_fastqc.html ${params.output}/dengue4/${x}/${x}_1_clean_fastqc.html
      mv ${params.output}/dengue4/${x}/${x}_1_fastqc.zip ${params.output}/dengue4/${x}/${x}_1_clean_fastqc.zip
      mv ${params.output}/dengue4/${x}/${x}_2_fastqc.html ${params.output}/dengue4/${x}/${x}_2_clean_fastqc.html
      mv ${params.output}/dengue4/${x}/${x}_2_fastqc.zip ${params.output}/dengue4/${x}/${x}_2_clean_fastqc.zip

   

   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
