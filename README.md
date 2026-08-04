# Daytona Dengue

<p align="center">
  <em>⚠️ For research use only. Results were obtained by procedures that were not CLIA validated.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Pipeline-Daytona%20Dengue-blue?style=plastic" />
  <img src="https://img.shields.io/badge/Nextflow-≥23.04-brightgreen?style=plastic&logo=nextflow" />
  <img src="https://img.shields.io/badge/Python-3.10+-yellow?style=plastic&logo=python" />
  <img src="https://img.shields.io/badge/License-MIT-red?style=plastic" />
</p>

## 🦟🧬 Overview

Daytona Dengue is Florida BPHL's Nextflow pipeline for Dengue virus (DENV) NGS data analysis. It processes paired-end Illumina reads through human read removal, quality control, adapter trimming, reference-based assembly, variant calling, serotype clade assignment and GenBank submission validation.

Serotype detection (DENV1–4) is performed automatically via Kraken2 and coverage-based screening and drives all downstream reference, primer and annotation selection. Nextclade provides fine-grained clade assignment using the [Hill et al. 2024 dengue lineage system (community/v-gen-lab datasets)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002834). VADR validates consensus sequences for GenBank submission.

### ⚙️ Dependencies

- **Nextflow** 23.04-26.x - [installation guide](https://github.com/nextflow-io/nextflow)
- **Apptainer/Singularity** - [installation guide](https://apptainer.org/docs/user/latest/)
- **Conda** - [installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- **SLURM** workload manager (required for HiPerGator; otherwise not required)

All bioinformatics tools run inside containers, no additional software installation is required.

### 💻 Resource Requirements

Daytona Dengue is designed to run on an HPC environment but can run locally with sufficient resources.

- **CPUs:** 24 recommended; minimum 8
- **RAM:** 50 GB recommended; minimum 16 GB
- **Disk:** ~2–3 GB per sample (input + output)

### 🛠️ Setup

#### 1. Clone this repository and enter the repository directory

```bash
$ git clone https://github.com/BPHL-Molecular/Daytona_dengue

$ cd Daytona_dengue/
```

#### 2. Create the conda environment

```bash
$ conda create -n daytona_dengue -c conda-forge python=3.10
```

#### 3. Configure params.yaml

Edit `params.yaml` and set the input and output paths for your run:

```yaml
input:  "/full/path/to/fastqs"
output: "/full/path/to/output"
```

Both `input` and `output` must be absolute paths with no trailing slash.

#### 4. Configure daytona_dengue.sh

> At Florida BPHL we use **Apptainer** on HiPerGator for containerization. `daytona_dengue.sh` is pre-configured for SLURM + Apptainer and is the recommended submission method for HiPerGator users.

Add your email address for job notifications and set `NXF_APPTAINER_CACHEDIR` to your Apptainer image cache directory:

```bash
export NXF_APPTAINER_CACHEDIR=/path/to/apptainer/cache
#SBATCH --mail-user=your@email.gov
```

### How to Run

Place paired FASTQ files in the directory specified by `params.input`. Both Illumina native (`SAMPLE_S1_L001_R1_001.fastq.gz`) and simplified (`SAMPLE_1.fastq.gz`) naming conventions are supported. If no matching FASTQ files are found, the pipeline exits immediately with an error.

### 🐊 HiPerGator Usage

```bash
sbatch daytona_dengue.sh
```

### ⚡ Local Usage

```bash
# Apptainer/Singularity
nextflow run daytona_dengue.nf -profile apptainer -params-file params.yaml
```

### Workflow Diagram

```mermaid
flowchart LR
    IN[Paired FASTQ] --> QC["Read QC and cleaning<br/>FastQC · Human Scrubber · Trimmomatic · BBTools"]
    QC --> SER["Serotype selection<br/>Kraken2 · BWA · Samtools coverage screen"]
    SER --> ASM["Assembly<br/>Samtools · iVar"]
    ASM --> VAL["Coverage QC and clade validation<br/>QC Gate · Nextclade · VADR"]

    QC --> REP[summary_report]
    SER --> REP
    ASM --> REP
    VAL --> REP

    VAL --> AQP[assemblies_qc_pass/]
    REP --> OUT[summary_report.txt]
    REP --> DASH[daytona_dengue_report.html]

    style SER fill:#9f9,stroke:#333,color:#000
    style REP fill:#f96,stroke:#333,stroke-width:2px,color:#000
    style OUT fill:#f96,stroke:#333,color:#000
    style DASH fill:#f96,stroke:#333,color:#000
    style AQP fill:#f96,stroke:#333,color:#000
```

> **QC GATE vs. assembly validation:** The **QC GATE** (`qc_flag`) is a minimum coverage (5%) and read-depth check (mean depth ≥ 30×) that confirms a sample's serotype classification is backed by enough on-target data, it is a serotype-classification QC, **not** a final assembly verdict, which is useful for surveillance purposes. Genome **assembly QC is performed by VADR**: a VADR **PASS** (`vadr_flag`) marks a submission-ready consensus, which is collected in `assemblies_qc_pass/`.

### 🧩 Modules

Daytona Dengue is made possible thanks to the following tools:

<small>

**Quality Control**: [FastQC](https://github.com/s-andrews/FastQC) 0.12.1 · [Trimmomatic](https://github.com/usadellab/Trimmomatic) 0.40 · [BBTools](https://github.com/bbushnell/BBTools) 39.84 · [MultiQC](https://github.com/MultiQC/MultiQC) 1.34

**Human Read Removal**: [NCBI SRA Human Scrubber](https://github.com/ncbi/sra-human-scrubber) 2.2.1

**Taxonomic Classification**: [Kraken2](https://github.com/DerrickWood/kraken2) 2.17.1

**Reference-Based Assembly**: [BWA](https://github.com/lh3/bwa) 0.7.19 · [Samtools](https://github.com/samtools/samtools) 1.23.1 · [iVar](https://github.com/andersen-lab/ivar) 1.4.4

**Clade Assignment**: [Nextclade](https://github.com/nextstrain/nextclade) 3.21.2 · [v-gen-lab dengue datasets](https://github.com/nextstrain/nextclade_data/tree/master/data/community/v-gen-lab/dengue) (Hill et al. 2024)

**Submission Validation**: [VADR](https://github.com/ncbi/vadr) 1.7

</small>

### 📁 Output

Per-sample results are written to `params.output/<sample_id>/`. A single summary file is written to `params.output/`:

```markdown
output/
├── <sample_id>/
│   ├── fastqc/
│   ├── trimmomatic/
│   ├── bbtools/
│   ├── bwa/
│   ├── samtools/
│   ├── ivar/
│   ├── vadr/
│   ├── nextclade/
│   └── multiqc/          # per-sample MultiQC (raw + clean FastQC)
├── assemblies_qc_pass/
├── daytona_dengue_report.html   # interactive run-level dashboard
└── summary_report.txt
```

| File | Samples | Key fields |
|------|---------|------------|
| `summary_report.txt` | All (including unclassified) | sample_id · serotype · nextclade_clade · kraken2_percent · reference · coverage stats · assembly stats · VADR_flag · QC_flag |
| `daytona_dengue_report.html` | All | Interactive dashboard: serotype/clade and coverage QC, assembly QC, raw and clean FastQC, software versions |
| `<sample_id>/multiqc/<sample_id>_multiqc_report.html` | Per sample | Raw + clean FastQC for that sample |

### 🤝 Contributing

We welcome contributions to make Daytona Dengue better! Feel free to open issues or submit pull requests to suggest additional features or enhancements.

### 📧 Contact

**Email**: [bphl-sebioinformatics@flhealth.gov](mailto:bphl-sebioinformatics@flhealth.gov)

### ⚖️ License

Daytona Dengue is licensed under the [MIT License](LICENSE).
