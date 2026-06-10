# Daytona Dengue

<p align="center">
  <em>⚠️ For research use only. Results were obtained by procedures that were not CLIA validated.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Pipeline-Daytona%20Dengue-blue?style=plastic" />
  <img src="https://img.shields.io/badge/Nextflow-≥23.04-brightgreen?style=plastic&logo=nextflow" />
  <img src="https://img.shields.io/badge/Python-3.10+-yellow?style=plastic&logo=python" />
  <img src="https://img.shields.io/badge/License-Apache%202.0-red?style=plastic" />
</p>

## 🦟🧬 Overview

Daytona Dengue is Florida BPHL's Nextflow pipeline for Dengue virus (DENV) NGS data analysis. It processes paired-end Illumina reads through human read removal, quality control, adapter trimming, reference-based assembly, variant calling, serotype clade assignment and GenBank submission validation.

Serotype detection (DENV1–4) is performed automatically via Kraken2 and coverage-based screening, and drives all downstream reference, primer, and annotation selection. Nextclade provides fine-grained clade assignment using the Hill et al. 2024 dengue lineage system (community/v-gen-lab datasets). VADR validates consensus sequences for GenBank submission.

### ⚙️ Dependencies

- **Nextflow** ≥ 23.04 — [installation guide](https://github.com/nextflow-io/nextflow)
- **Apptainer/Singularity** — [installation guide](https://apptainer.org/docs/user/latest/)
- **SLURM** workload manager (required for HiPerGator; optional otherwise)

All bioinformatics tools run inside containers — no additional software installation is required.

### 💻 Resource Requirements

Daytona Dengue is designed to run on an HPC environment but can run locally with sufficient resources.

- **CPUs:** 40 recommended (HPC); minimum 8
- **RAM:** 200 GB recommended (HPC); 32 GB minimum (Kraken2 requires ~32 GB alone)
- **Disk:** ~5 GB per sample (input + output); ~8–32 GB for the Kraken2 viral database

**Estimated runtime** (18 samples, 40 CPUs, HPC): ~2–3 hours, dominated by Kraken2 and BWA alignment.

### 🛠️ Setup

#### 1. Configure params.yaml

Edit `params.yaml` and set the input and output paths for your run:

```yaml
input:  "/full/path/to/fastqs"
output: "/full/path/to/output"
```

Both `input` and `output` must be absolute paths with no trailing slash.

> **HiPerGator users:** only `input` and `output` need to be set. All other paths are pre-configured in `nextflow.config`.
> **Non-HiPerGator users:** set `params.kraken_db` in `nextflow.config` (or add it to `params.yaml`) to point to your local Kraken2 viral database directory.

#### 2. Configure daytona_dengue.sh

Set `NXF_APPTAINER_CACHEDIR` to your Apptainer image cache directory and add your email address for job notifications:

```bash
export NXF_APPTAINER_CACHEDIR=/path/to/apptainer/cache
#SBATCH --mail-user=your@email.gov
```

### How to Run

Place paired FASTQ files in the directory specified by `params.input`. Both Illumina native (`SAMPLE_S1_L001_R1_001.fastq.gz`) and simplified (`SAMPLE_1.fastq.gz`) naming conventions are supported.

### 🐊 HiPerGator Usage

```bash
sbatch daytona_dengue.sh
```

### ⚡ Local Usage

```bash
# Apptainer/Singularity
nextflow run daytona_dengue.nf -profile apptainer -params-file params.yaml

# Docker
nextflow run daytona_dengue.nf -profile docker -params-file params.yaml
```

### Workflow Diagram

```mermaid
flowchart LR
    A[Paired FASTQ Input] --> B[FastQC]
    A --> HS[Human Scrubber]
    HS --> C[Trimmomatic]
    C --> D[BBTools\nadapters]
    D --> E[BBTools\nPhiX]
    E --> F[FastQC clean]
    B --> MQ[MultiQC]
    F --> MQ

    E --> K[Kraken2]
    E --> BWA[BWA\n4× DENV refs]
    BWA --> SS[Samtools screen\ncoverage]
    SS --> SD[Serotype Detect]

    SD --> |DENV1-4| SB[Samtools BAM\nwinning ref]
    SB --> IT[iVar trim]
    IT --> SC[Samtools coverage]
    IT --> SM[Samtools mpileup]
    SM --> IV[iVar variants]
    SM --> IC[iVar consensus]
    IC --> QG[QC Gate]

    QG --> |PASS| VD[VADR]
    QG --> |PASS| NC[Nextclade\nper-serotype dataset]

    K --> SR[summary_report]
    SD --> SR
    SC --> SR
    IC --> SR
    NC --> SR
    VD --> SR
    QG --> SR
    SR --> OUT[summary_report.txt]

    style SD fill:#fef,stroke:#333,color:#000
    style NC fill:#9cf,stroke:#333,color:#000
    style VD fill:#9cf,stroke:#333,color:#000
    style SR fill:#f96,stroke:#333,stroke-width:2px,color:#000
    style OUT fill:#f96,stroke:#333,stroke-width:3px,color:#000
```

### 🧩 Modules

Daytona Dengue is made possible thanks to the following tools:

<small>

**Quality Control** — [FastQC](https://github.com/s-andrews/FastQC) 0.12.1 · [Trimmomatic](https://github.com/usadellab/Trimmomatic) 0.40 · [BBTools](https://github.com/bbushnell/BBTools) 39.84 · [MultiQC](https://github.com/MultiQC/MultiQC) 1.34

**Human Read Removal** — [NCBI SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) 2.2.1

**Taxonomic Classification** — [Kraken2](https://github.com/DerrickWood/kraken2) 2.17.1

**Reference-Based Assembly** — [BWA](https://github.com/lh3/bwa) 0.7.19 · [Samtools](https://github.com/samtools/samtools) 1.23.1 · [iVar](https://github.com/andersen-lab/ivar) 1.4.4

**Clade Assignment** — [Nextclade](https://github.com/nextstrain/nextclade) 3.21.2 · [v-gen-lab dengue datasets](https://github.com/nextstrain/nextclade_data/tree/master/data/community/v-gen-lab/dengue) (Hill et al. 2024)

**Submission Validation** — [VADR](https://github.com/ncbi/vadr) 1.7

</small>

### 📁 Output

Per-sample results are written to `params.output/<sample_id>/`. A single summary file is written to `params.output/`:

```
output/
├── <sample_id>/
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bbtools/
│   ├── bwa/
│   ├── samtools/
│   ├── ivar/
│   ├── vadr/
│   └── nextclade/
├── multiqc/
└── summary_report.txt
```

| File | Samples | Key fields |
|------|---------|------------|
| `summary_report.txt` | All (including unclassified) | sample_id · serotype · nextclade_clade · nextclade_qc_overall · kraken2_percent · reference · coverage stats · assembly stats · VADR_flag · QC_flag |

### 📁 Directory Structure

```text
daytona_dengue/
├── daytona_dengue.nf       # Entry workflow
├── daytona_dengue.sh       # SLURM submission script
├── nextflow.config         # Container, resource, and profile configuration
├── params.yaml             # Runtime parameters (input/output paths)
├── modules/                # One .nf file per tool or tool group
├── bin/                    # Python helper scripts (QC gate, summary report)
└── assets/
    ├── reference/          # Reference FASTA + BWA indexes (DENV1–4)
    ├── primers/            # Primer BED files (DENV1–4)
    ├── annotations/        # Per-serotype GFF files (for iVar variants)
    ├── nextclade/          # Cached Nextclade datasets (auto-downloaded)
    └── vadr/               # Cached VADR models (auto-downloaded)
```

### 🤝 Contributing

We welcome contributions to make Daytona Dengue better! Feel free to open issues or submit pull requests to suggest additional features or enhancements.

### 📧 Contact

**Email**: [bphl-sebioinformatics@flhealth.gov](mailto:bphl-sebioinformatics@flhealth.gov)

### ⚖️ License

Daytona Dengue is licensed under the [Apache License 2.0](LICENSE).
