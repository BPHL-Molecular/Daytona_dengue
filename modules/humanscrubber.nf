process humanscrubber {
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

      # Run sra-human-scrubber to remove human reads
      gzip -d ${params.output}/dengue1/${x}/${x}_1.fastq.gz
      gzip -d ${params.output}/dengue1/${x}/${x}_2.fastq.gz

      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue1/${x}/${x}_1.fastq -o ${params.output}/dengue1/${x}/${x}_1_humanclean.fastq
      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue1/${x}/${x}_2.fastq -o ${params.output}/dengue1/${x}/${x}_2_humanclean.fastq

      gzip ${params.output}/dengue1/${x}/${x}_1_humanclean.fastq
      gzip ${params.output}/dengue1/${x}/${x}_2_humanclean.fastq

      
   elif [[ ${x} =~ SER2_ ]];then

      # Run sra-human-scrubber to remove human reads
      gzip -d ${params.output}/dengue2/${x}/${x}_1.fastq.gz
      gzip -d ${params.output}/dengue2/${x}/${x}_2.fastq.gz

     /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue2/${x}/${x}_1.fastq -o ${params.output}/dengue2/${x}/${x}_1_humanclean.fastq
     /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue2/${x}/${x}_2.fastq -o ${params.output}/dengue2/${x}/${x}_2_humanclean.fastq

      gzip ${params.output}/dengue2/${x}/${x}_1_humanclean.fastq
      gzip ${params.output}/dengue2/${x}/${x}_2_humanclean.fastq

      
   elif [[ ${x} =~ SER3_ ]];then

      # Run sra-human-scrubber to remove human reads
      gzip -d ${params.output}/dengue3/${x}/${x}_1.fastq.gz
      gzip -d ${params.output}/dengue3/${x}/${x}_2.fastq.gz

      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue3/${x}/${x}_1.fastq -o ${params.output}/dengue3/${x}/${x}_1_humanclean.fastq
      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue3/${x}/${x}_2.fastq -o ${params.output}/dengue3/${x}/${x}_2_humanclean.fastq

      gzip ${params.output}/dengue3/${x}/${x}_1_humanclean.fastq
      gzip ${params.output}/dengue3/${x}/${x}_2_humanclean.fastq
   

      
   elif [[ ${x} =~ SER4_ ]];then

      # Run sra-human-scrubber to remove human reads
      gzip -d ${params.output}/dengue4/${x}/${x}_1.fastq.gz
      gzip -d ${params.output}/dengue4/${x}/${x}_2.fastq.gz

      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue4/${x}/${x}_1.fastq -o ${params.output}/dengue4/${x}/${x}_1_humanclean.fastq
      /opt/scrubber/scripts/scrub.sh -r -i ${params.output}/dengue4/${x}/${x}_2.fastq -o ${params.output}/dengue4/${x}/${x}_2_humanclean.fastq

      gzip ${params.output}/dengue4/${x}/${x}_1_humanclean.fastq
      gzip ${params.output}/dengue4/${x}/${x}_2_humanclean.fastq

   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
