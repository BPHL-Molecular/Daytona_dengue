# Daytona_dengue
A nextflow pipeline for Dengue NGS data analysis. 

## What to do
The pipeline can analyze NGS data of Dengue virus. The sample's serotype, sequencing quality, mapping/alignment with reference, coverage, SNPs, annotation, etc can be outputted.  

## Prerequisites
Nextflow should be installed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed.

Singularity/Apptainer is needed if you select singularity to run containers in the pipeline. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.

Docker ( https://www.docker.com/ ) is needed if you select docker to run containers in the pipeline.

SLURM is needed if you plan use slurm to load the pipeline to run at high performance computing cluster. 

## Before running
If the referene genomes donot have index, using below command to generate their index files and to put them in the same directory before running the pipeline.
```bash
singularity exec docker://staphb/bwa:0.7.17 bwa index <full path to your genome fasta file> 
```
         
## How to run
### By default: 
By default, the pipeline uses singularity to run containers and is wrapped by SLURM. To run the pipeline by default, you should follow the steps below.
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz".
2. open file "parames.yaml", set the parameters. 
3. get into the top of the pipeline directory, then run 
```bash
sbatch ./daytona_dengue.sh
```
Note: the sbatch parameters setting is based on our cluser HiPerGator at University of Florida campus. You maybe need change them according to your cluster's configuration.
### By Singularity
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz".
2. open file "parames.yaml", set the parameters. 
3. get into the top of the pipeline directory, then run 
```bash
bash ./kraken2_viral.sh
``` 
```bash
nextflow run daytona_dengue.nf -params-file params.yaml -c ./configs/singularity.config
``` 
### By Docker: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz".
2. open file "parames.yaml", set the parameters. 
3. get into the top of the pipeline directory, then run 
```bash
bash ./kraken2_viral.sh
```
```bash
nextflow run daytona_dengue.nf -params-file params.yaml -c ./configs/docker.config
```
## Results
All results can be found in the directory /output.

## Note
1. In the output file Serotypes.txt, a sample is considered as unserotypeed if its confident rate under 50%. This sample's exact taxon can be found in /output/kraken_out_broad.
2. renamefile.sh can be used to change user's file name
