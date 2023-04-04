process trimmomatic {
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
      #Run trimmomatic
      trimmomatic PE -phred33 -trimlog ${params.output}/dengue1/${x}/${x}.log ${params.output}/dengue1/${x}/${x}_1_humanclean.fastq.gz ${params.output}/dengue1/${x}/${x}_2_humanclean.fastq.gz ${params.output}/dengue1/${x}/${x}_trim_1.fastq.gz ${params.output}/dengue1/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/dengue1/${x}/${x}_trim_2.fastq.gz ${params.output}/dengue1/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/dengue1/${x}/${x}_trimstats.txt

      rm ${params.output}/dengue1/${x}/${x}_unpaired_trim_*.fastq.gz
      rm ${params.output}/dengue1/${x}/${x}_1.fastq ${params.output}/dengue1/${x}/${x}_2.fastq

      
   elif [[ ${x} =~ SER2_ ]];then
      #Run trimmomatic
      trimmomatic PE -phred33 -trimlog ${params.output}/dengue2/${x}/${x}.log ${params.output}/dengue2/${x}/${x}_1_humanclean.fastq.gz ${params.output}/dengue2/${x}/${x}_2_humanclean.fastq.gz ${params.output}/dengue2/${x}/${x}_trim_1.fastq.gz ${params.output}/dengue2/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/dengue2/${x}/${x}_trim_2.fastq.gz ${params.output}/dengue2/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/dengue2/${x}/${x}_trimstats.txt

      rm ${params.output}/dengue2/${x}/${x}_unpaired_trim_*.fastq.gz
      rm ${params.output}/dengue2/${x}/${x}_1.fastq ${params.output}/dengue2/${x}/${x}_2.fastq


      
   elif [[ ${x} =~ SER3_ ]];then
      #Run trimmomatic
      trimmomatic PE -phred33 -trimlog ${params.output}/dengue3/${x}/${x}.log ${params.output}/dengue3/${x}/${x}_1_humanclean.fastq.gz ${params.output}/dengue3/${x}/${x}_2_humanclean.fastq.gz ${params.output}/dengue3/${x}/${x}_trim_1.fastq.gz ${params.output}/dengue3/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/dengue3/${x}/${x}_trim_2.fastq.gz ${params.output}/dengue3/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/dengue3/${x}/${x}_trimstats.txt

      rm ${params.output}/dengue3/${x}/${x}_unpaired_trim_*.fastq.gz
      rm ${params.output}/dengue3/${x}/${x}_1.fastq ${params.output}/dengue3/${x}/${x}_2.fastq

   

      
   elif [[ ${x} =~ SER4_ ]];then
      #Run trimmomatic
      trimmomatic PE -phred33 -trimlog ${params.output}/dengue4/${x}/${x}.log ${params.output}/dengue4/${x}/${x}_1_humanclean.fastq.gz ${params.output}/dengue4/${x}/${x}_2_humanclean.fastq.gz ${params.output}/dengue4/${x}/${x}_trim_1.fastq.gz ${params.output}/dengue4/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/dengue4/${x}/${x}_trim_2.fastq.gz ${params.output}/dengue4/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/dengue4/${x}/${x}_trimstats.txt

      rm ${params.output}/dengue4/${x}/${x}_unpaired_trim_*.fastq.gz
      rm ${params.output}/dengue4/${x}/${x}_1.fastq ${params.output}/dengue4/${x}/${x}_2.fastq

   

   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
