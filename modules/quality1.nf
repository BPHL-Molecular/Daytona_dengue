process quality1 {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${params.output}/dengue1/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #echo ${params.input}/${x}_1.fastq.gz >> xfile.txt
   
   mkdir -p ${params.output}/dengue1
   mkdir -p ${params.output}/dengue1/assemblies
   mkdir -p ${params.output}/dengue1/variants
   mkdir -p ${params.output}/dengue1/vadr_error_reports
   mkdir -p ${params.output}/dengue1/${x}
   cp ${params.input}/dengue1/${x}_*.fastq.gz ${params.output}/dengue1/${x}
   
   #Run fastqc on original reads
   singularity exec --cleanenv docker://staphb/fastqc:0.11.9 fastqc ${params.output}/dengue1/${x}/${x}_1.fastq.gz ${params.output}/dengue1/${x}/${x}_2.fastq.gz

   mv ${params.output}/dengue1/${x}/${x}_1_fastqc.html ${params.output}/dengue1/${x}/${x}_1_original_fastqc.html
   mv ${params.output}/dengue1/${x}/${x}_1_fastqc.zip ${params.output}/dengue1/${x}/${x}_1_original_fastqc.zip
   mv ${params.output}/dengue1/${x}/${x}_2_fastqc.html ${params.output}/dengue1/${x}/${x}_2_original_fastqc.html
   mv ${params.output}/dengue1/${x}/${x}_2_fastqc.zip ${params.output}/dengue1/${x}/${x}_2_original_fastqc.zip
   
   # Run sra-human-scrubber to remove human reads
   gzip -d ${params.output}/dengue1/${x}/${x}_1.fastq.gz
   gzip -d ${params.output}/dengue1/${x}/${x}_2.fastq.gz

   singularity exec -B ${params.output}/dengue1/${x}:/data docker://ncbi/sra-human-scrubber:1.1.2021-05-05 /opt/scrubber/scripts/scrub.sh -r -i /data/${x}_1.fastq -o /data/${x}_1_humanclean.fastq
   singularity exec -B ${params.output}/dengue1/${x}:/data docker://ncbi/sra-human-scrubber:1.1.2021-05-05 /opt/scrubber/scripts/scrub.sh -r -i /data/${x}_2.fastq -o /data/${x}_2_humanclean.fastq

   gzip ${params.output}/dengue1/${x}/${x}_1_humanclean.fastq
   gzip ${params.output}/dengue1/${x}/${x}_2_humanclean.fastq
   
   #Run trimmomatic
   #singularity exec --cleanenv /apps/staphb-toolkit/containers/trimmomatic_0.39.sif trimmomatic PE -phred33 -trimlog ${params.output}/${x}/${x}.log ${params.output}/${x}/${x}_1_humanclean.fastq.gz ${params.output}/${x}/${x}_2_humanclean.fastq.gz ${params.output}/${x}/${x}_trim_1.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/${x}/${x}_trim_2.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/${x}/${x}_trimstats.txt
   singularity exec --cleanenv docker://staphb/trimmomatic:0.39 trimmomatic PE -phred33 -trimlog ${params.output}/dengue1/${x}/${x}.log ${params.output}/dengue1/${x}/${x}_1_humanclean.fastq.gz ${params.output}/dengue1/${x}/${x}_2_humanclean.fastq.gz ${params.output}/dengue1/${x}/${x}_trim_1.fastq.gz ${params.output}/dengue1/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/dengue1/${x}/${x}_trim_2.fastq.gz ${params.output}/dengue1/${x}/${x}_unpaired_trim_2.fastq.gz SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20 > ${params.output}/dengue1/${x}/${x}_trimstats.txt

   rm ${params.output}/dengue1/${x}/${x}_unpaired_trim_*.fastq.gz
   rm ${params.output}/dengue1/${x}/${x}_1.fastq ${params.output}/dengue1/${x}/${x}_2.fastq

   #Run bbduk to remove Illumina adapter sequences and any PhiX contamination  
   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/dengue1/${x}/${x}_trim_1.fastq.gz in2=${params.output}/dengue1/${x}/${x}_trim_2.fastq.gz out1=${params.output}/dengue1/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/dengue1/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/dengue1/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/dengue1/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/dengue1/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/dengue1/${x}/${x}_1.fq.gz out2=${params.output}/dengue1/${x}/${x}_2.fq.gz outm=${params.output}/dengue1/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/dengue1/${x}/${x}_phixstats.txt

   rm ${params.output}/dengue1/${x}/${x}_trim*.fastq.gz
   rm ${params.output}/dengue1/${x}/${x}*rmadpt.fq.gz
   
   #Run fastqc on clean forward and reverse reads
   singularity exec --cleanenv docker://staphb/fastqc:0.11.9 fastqc ${params.output}/dengue1/${x}/${x}_1.fq.gz ${params.output}/dengue1/${x}/${x}_2.fq.gz

   #Rename fastqc output files
   mv ${params.output}/dengue1/${x}/${x}_1_fastqc.html ${params.output}/dengue1/${x}/${x}_1_clean_fastqc.html
   mv ${params.output}/dengue1/${x}/${x}_1_fastqc.zip ${params.output}/dengue1/${x}/${x}_1_clean_fastqc.zip
   mv ${params.output}/dengue1/${x}/${x}_2_fastqc.html ${params.output}/dengue1/${x}/${x}_2_clean_fastqc.html
   mv ${params.output}/dengue1/${x}/${x}_2_fastqc.zip ${params.output}/dengue1/${x}/${x}_2_clean_fastqc.zip
   
   #Run multiqc
   singularity exec --cleanenv docker://staphb/multiqc:1.8 multiqc ${params.output}/dengue1/${x}/${x}_*_fastqc.zip -o ${params.output}/dengue1/${x}

   #Map reads to reference
   mkdir ${params.output}/dengue1/${x}/alignment
   """
}
