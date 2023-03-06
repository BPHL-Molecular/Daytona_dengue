# Daytona_dengue
A nextflow pipeline for Dengue NGS data analysis. 

## What to do
The pipeline can analyze NGS data of Dengue virus. The sample's serotype, sequencing quality, mapping/alignment with reference, coverage, SNPs, annotation, etc can be outputted.  

## Prerequisites
Nextflow should be installed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed.

Singularity/Apptainer is also needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.

## Before running
If the referene genomes donot have index, using below command to generate their index files and to put them in the same directory before running the pipeline.
         ```singularity exec docker://staphb/bwa:0.7.17 bwa index <full path to your genome fasta file> ```
         
## How to run
### 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz".
2. open file "parames.yaml", set the parameters. 
3. get into the directory of the pipeline, run "sbatch ./daytona_dengue.sh"

## Results
All results can be found in the directory /output.

