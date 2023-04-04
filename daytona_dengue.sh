#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=Daytona_dengue
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb
#SBATCH --output=daytona_dengue.%j.out
#SBATCH --error=daytona_dengue.%j.err
#SBATCH --time=3-00

#module load singularity
module load apptainer

#identify serotype of each sample
bash ./kraken2_viral.sh

#nextflow run daytona_dengue.nf -params-file params.yaml -c ./configs/singularity.config
#nextflow run daytona_dengue.nf -params-file params.yaml -c ./configs/docker.config
nextflow run daytona_dengue.nf -params-file params.yaml


sort ./output/dengue*/*/report.txt | uniq > ./output/sum_report.txt
sed -i '/sampleID\treference/d' ./output/sum_report.txt
sed -i '1i sampleID\treference\tstart\tend\tnum_raw_reads\tnum_clean_reads\tnum_mapped_reads\tpercent_mapped_clean_reads\tcov_bases_mapped\tpercent_genome_cov_map\tmean_depth\tmean_base_qual\tmean_map_qual\tassembly_length\tnumN\tpercent_ref_genome_cov\tVADR_flag\tQC_flag' ./output/sum_report.txt

mv ./Serotypes.txt ./output/
mv ./kraken_out_broad ./output/

#cat ./output/assemblies/*.fa > ./output/assemblies.fasta
#singularity exec /apps/staphb-toolkit/containers/nextclade_2021-03-15.sif nextclade --input-fasta ./output/assemblies.fasta --output-csv ./output/nextclade_report.csv
