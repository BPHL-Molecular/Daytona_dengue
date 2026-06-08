# Daytona Dengue

A Nextflow DSL2 pipeline for Dengue virus NGS data analysis developed at Florida's Bureau of Public Health Laboratories (BPHL).

The pipeline processes paired-end Illumina reads through quality control, human read removal, adapter trimming, reference-based assembly, variant calling, and sequence annotation. Serotype detection (DENV1–4) is performed inside the pipeline via Kraken2, driving all downstream reference and primer selection automatically. Nextclade provides clade assignment; VADR validates sequences for GenBank submission.

## Directory Structure

```
daytona_dengue/
├── daytona_dengue.nf       # Entry workflow
├── daytona_dengue.sh       # SLURM submission script
├── nextflow.config         # All configuration — params, containers, resources, profiles
├── params.yaml             # Runtime parameters
├── CHANGELOG.md
├── README.md
├── modules/                # One .nf file per tool or tool group
├── bin/                    # Python helper scripts
├── assets/
│   └── annotations/        # Per-serotype GFF annotation files (for iVar)
├── reference/              # Reference FASTA files + BWA indexes (DENV1–4)
├── primers/                # Primer BED files (DENV1–4)
└── fastqs/                 # Place input FASTQ files here
```

## Prerequisites

- [Nextflow](https://github.com/nextflow-io/nextflow) ≥ 23.04
- [Apptainer](https://apptainer.org/) (recommended for HPC) or Docker
- SLURM (for HPC runs)
- A Kraken2 viral database (path set in `params.yaml`)

## How to Run

### HPC (SLURM + Apptainer)

1. Place paired FASTQ files in `fastqs/`. Files must match either `*_{1,2}.fastq.gz` or `*_R{1,2}_*.fastq.gz` naming.
2. Edit `params.yaml` — set `input`, `output`, `reference`, `primer`, `annotations`, and `kraken_db` paths.
3. Edit `daytona_dengue.sh` — set `--mail-user` and `conda activate <your-env>`.
4. Submit:
```bash
sbatch daytona_dengue.sh
```

### Local (Docker)

```bash
nextflow run daytona_dengue.nf -profile docker -params-file params.yaml
```

### Local (Singularity/Apptainer)

```bash
nextflow run daytona_dengue.nf -profile apptainer -params-file params.yaml
```

## Results

Per-sample results are written to subdirectories under the `output` path defined in `params.yaml`:

```
output/
├── <sample_id>/
│   ├── kraken2/
│   ├── fastqc/
│   ├── humanscrubber/
│   ├── trimmomatic/
│   ├── bbtools/
│   ├── bwa/
│   ├── samtools/
│   ├── ivar/
│   ├── stats/
│   ├── vadr/
│   └── nextclade/
├── multiqc/
└── summary_report.tsv      # Final aggregated report
```

`summary_report.tsv` contains per-sample serotype, QC metrics, Nextclade clade assignment, and VADR submission result.

## Reference Data

| File type | Location |
|---|---|
| Reference FASTA + BWA indexes | `reference/` |
| Primer BED files | `primers/` |
| GFF annotation files (iVar) | `assets/annotations/` |
| Kraken2 database | External path — set in `params.yaml` |
| Nextclade dataset | Auto-downloaded and cached (see `nextclade_cache_dir` in `params.yaml`) |
