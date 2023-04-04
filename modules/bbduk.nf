process bbduk {
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
      #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
      bbduk.sh in1=${params.output}/dengue1/${x}/${x}_trim_1.fastq.gz in2=${params.output}/dengue1/${x}/${x}_trim_2.fastq.gz out1=${params.output}/dengue1/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/dengue1/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/dengue1/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
      bbduk.sh in1=${params.output}/dengue1/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/dengue1/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/dengue1/${x}/${x}_1.fq.gz out2=${params.output}/dengue1/${x}/${x}_2.fq.gz outm=${params.output}/dengue1/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/dengue1/${x}/${x}_phixstats.txt

      rm ${params.output}/dengue1/${x}/${x}_trim*.fastq.gz
      rm ${params.output}/dengue1/${x}/${x}*rmadpt.fq.gz

      
   elif [[ ${x} =~ SER2_ ]];then
      #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
      bbduk.sh in1=${params.output}/dengue2/${x}/${x}_trim_1.fastq.gz in2=${params.output}/dengue2/${x}/${x}_trim_2.fastq.gz out1=${params.output}/dengue2/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/dengue2/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/dengue2/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
      bbduk.sh in1=${params.output}/dengue2/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/dengue2/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/dengue2/${x}/${x}_1.fq.gz out2=${params.output}/dengue2/${x}/${x}_2.fq.gz outm=${params.output}/dengue2/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/dengue2/${x}/${x}_phixstats.txt

      rm ${params.output}/dengue2/${x}/${x}_trim*.fastq.gz
      rm ${params.output}/dengue2/${x}/${x}*rmadpt.fq.gz


   elif [[ ${x} =~ SER3_ ]];then
      #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
      bbduk.sh in1=${params.output}/dengue3/${x}/${x}_trim_1.fastq.gz in2=${params.output}/dengue3/${x}/${x}_trim_2.fastq.gz out1=${params.output}/dengue3/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/dengue3/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/dengue3/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
      bbduk.sh in1=${params.output}/dengue3/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/dengue3/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/dengue3/${x}/${x}_1.fq.gz out2=${params.output}/dengue3/${x}/${x}_2.fq.gz outm=${params.output}/dengue3/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/dengue3/${x}/${x}_phixstats.txt

      rm ${params.output}/dengue3/${x}/${x}_trim*.fastq.gz
      rm ${params.output}/dengue3/${x}/${x}*rmadpt.fq.gz

  
   elif [[ ${x} =~ SER4_ ]];then
      #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
      bbduk.sh in1=${params.output}/dengue4/${x}/${x}_trim_1.fastq.gz in2=${params.output}/dengue4/${x}/${x}_trim_2.fastq.gz out1=${params.output}/dengue4/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/dengue4/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/dengue4/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
      bbduk.sh in1=${params.output}/dengue4/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/dengue4/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/dengue4/${x}/${x}_1.fq.gz out2=${params.output}/dengue4/${x}/${x}_2.fq.gz outm=${params.output}/dengue4/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/dengue4/${x}/${x}_phixstats.txt

      rm ${params.output}/dengue4/${x}/${x}_trim*.fastq.gz
      rm ${params.output}/dengue4/${x}/${x}*rmadpt.fq.gz


   else
      echo "No serotyped sequence is in fastqs folder"
   fi

   """
}
