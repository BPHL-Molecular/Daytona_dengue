# Changelog

All notable changes to the Daytona Dengue pipeline will be documented in this file.

---

## [Unreleased] — Post-refactor improvements

### Changed

- `modules/kraken2.nf` — reverted to bundled viral container (`staphb/kraken2:2.17.1-viral-20251015`); removed `--confidence` filter
- `nextflow.config` — reduced Kraken2 memory from 50 GB to 8 GB
- `modules/nextclade.nf` — switched to per-serotype community datasets (`community/v-gen-lab/dengue/denv1–4`)
- `bin/summary_report.py` — renamed `VADR_flag`/`QC_flag` → `vadr_flag`/`qc_flag`; removed `nextclade_qc_overall` column; improved unclassified/low-coverage QC message; VADR flag derived from QC flag for samples that never reached VADR
- `bin/qc_gate.py` — updated QC coverage threshold; depth fail message now shows actual observed depth
- `README.md` — corrected resource requirements; added conda setup step; added Nextflow version constraint note (≥23.04, <26.0); removed estimated runtime
- `daytona_dengue.sh` — updated SLURM resource allocation

---

## [Unreleased] — SOP Format Refactor

### Added

- `bin/` directory for Python helper scripts
- `assets/annotations/` directory for GFF annotation files
- `modules/kraken2.nf` — Kraken2 serotyping as a Nextflow process
- `modules/vadr.nf` — VADR annotation as a dedicated process (QC-pass samples only)
- `modules/nextclade_download.nf` — Nextclade dataset download with storeDir caching
- `modules/nextclade.nf` — Nextclade clade assignment on all consensus sequences
- `modules/summary_report.nf` — barrier-synchronized final report generation
- `bin/parse_kraken2.py` — parses Kraken2 report to extract serotype
- `bin/pystats.py` — QC metrics script (extracted from pystats.nf)
- `bin/summary_report.py` — aggregates per-sample stats into final TSV
- `CHANGELOG.md`

### Changed

- `nextflow.config` — full rewrite; single config replaces `configs/` directory; all processes get `container`, `cpus`, `memory`; parameterized tool settings; profiles for `standard`, `docker`, `singularity`, `apptainer`
- `params.yaml` — updated with `annotations`, `kraken_db`, `nextclade_cache_dir` keys; tool params as commented-out overrides
- `daytona_dengue.sh` — updated to BPHL guide format; post-processing removed; timestamp rename block added
- `daytona_dengue.nf` — full rewrite; `meta`-driven channel design; `meta.serotype` replaces all `if/elif SER1–4` blocks
- `modules/fastqc.nf` — modernized; absorbs `fastqc_clean.nf`
- `modules/humanscrubber.nf` — modernized; `tag`, `publishDir`, named emits
- `modules/trimmomatic.nf` — modernized; params from config; `task.cpus`
- `modules/bbtools.nf` — renamed from `bbduk.nf`; split into `bbtools_adapters` + `bbtools_phix`
- `modules/bwa.nf` — renamed from `frag_bwa.nf`; reference passed as channel input
- `modules/multiqc.nf` — modernized
- `modules/samtools.nf` — merged from `frag_samtools.nf` + `primer_trim_samtools.nf`
- `modules/ivar.nf` — merged from `primer_trim_ivar.nf` + `assembly.nf`
- `modules/pystats.nf` — thin wrapper; all logic moved to `bin/pystats.py`
- GFF files moved from `reference/` to `assets/annotations/`

### Removed

- `configs/` directory — replaced by single `nextflow.config`
- `modules/fastqc_clean.nf` — merged into `modules/fastqc.nf`
- `modules/frag_bwa.nf` — renamed to `modules/bwa.nf`
- `modules/frag_samtools.nf` — merged into `modules/samtools.nf`
- `modules/primer_trim_samtools.nf` — merged into `modules/samtools.nf`
- `modules/primer_trim_ivar.nf` — merged into `modules/ivar.nf`
- `modules/assembly.nf` — merged into `modules/ivar.nf`
- `kraken2_viral.sh` — replaced by `modules/kraken2.nf`
- `renamefile.sh` — replaced by `id.replaceAll()` in input channel
- `table.py` — replaced by `bin/summary_report.py`
- `rename_aqp.py` — replaced by meta-driven design
