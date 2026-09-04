# Changelog

All notable changes to the Daytona Dengue pipeline will be documented in this file.

---

## [v1.0.0] - SOP Format Refactor, Custom MultiQC dashboard and Nextflow 26 support

### Added

- `assets/multiqc_config.yaml` and `assets/daytona_dengue_report.css` - branded, interactive
  run-level MultiQC dashboard (`daytona_dengue_report.html`), replacing the generic aggregate
  report previously published to `output/all_multiqc/`
- `bin/summary_report.py` - `_mqc_preamble`/`_write_mqc`/`emit_daytona_mqc_tables` write
  `daytona_dengue_serotype_mqc.tsv` and `daytona_dengue_assembly_mqc.tsv`, rendered as
  sortable, color-coded tables (`qc_flag`, `vadr_flag`, `serotype` verdicts) in the dashboard
- `modules/multiqc.nf` - `multiqc` process now stages the config/CSS/tables, writes an inline
  Software Versions table from the container tags pinned in `nextflow.config`, and cleans up
  the `_mqc.tsv` files after each run so a resumed run doesn't re-ingest stale tables
- `manifest { nextflowVersion = '>=23.04' }` in `nextflow.config`
- `modules/fastqc.nf` - `fastqc` process symlinks its inputs to `<sample>_R{1,2}_raw.fastq.gz`
  before running, so raw and clean reads collapse onto one General Statistics row per sample
  instead of two. MultiQC keys FastQC rows off the `Filename` recorded inside `fastqc_data.txt`,
  which follows the input name, so renaming the output zip has no effect

### Changed

- All 23 `publishDir` directives across `modules/*.nf` rewritten from bare-string to closure
  form (`publishDir { "..." }, mode: 'copy'`), required by the v2 strict script parser that
  Nextflow 26.04 makes the default. The bare-string form evaluates the path at process-definition
  time, before `meta` is in scope, raising `No such variable: meta` at module load
- `daytona_dengue.nf` - `ch_barrier` gets `.ifEmpty(true)`, so a run where no sample reaches
  `vadr`/`nextclade` (e.g. every sample fails serotyping) still triggers `summary_report`
  instead of stalling with no error
- `daytona_dengue.sh` - loads the default `nextflow` module instead of pinning `nextflow/25.10.4`
- `assets/multiqc_config.yaml` - FastQC now runs twice (raw and clean reads), narrowing General
  Statistics to six ordered metrics (Seqs, Median len, GC for each) and dropping the BBTools,
  Kraken, Samtools and Trimmomatic sections and columns via `exclude_modules`; DENV1-4 serotype
  calls colored green in `table_cond_formatting_rules`, matching the PASS coloring
- `bin/summary_report.py` - `'Serotype ID and Coverage QC'` renamed to
  `'Serotype/Clade and Coverage QC'` and gains `nextclade_clade`, dropping `kraken2_percent`;
  `'Assembly and Clade QC'` renamed to `'Assembly QC'` and loses `nextclade_clade`
- `README.md` - Nextflow support range updated to 23.04-26.x; output section reflects
  `daytona_dengue_report.html` replacing `all_multiqc/`, and the dashboard description reflects
  the renamed sections and dual FastQC columns
- `bin/summary_report.py` and `bin/qc_gate.py` - `percent_genome_cov_map` renamed to
  `percent_genome_cov_aligned` (breadth of coverage from the mapped BAM) and
  `percent_ref_genome_cov` renamed to `percent_genome_cov_assembled` (completeness of the
  final iVar consensus, the value `qc_flag` is actually thresholded on); the
  Serotype/Clade and Coverage QC table now shows `percent_genome_cov_assembled` next to
  `qc_flag` instead of `percent_genome_cov_aligned`, so the displayed number matches the
  verdict it drives
  
### Added

- Per-sample MultiQC (`multiqc_sample`) over each sample's raw + clean FastQC, published to
  `output/<sample_id>/multiqc/`

### Changed

- Aggregate MultiQC output dir renamed to `output/all_multiqc/`; aggregate now runs with
  `--ignore "*/multiqc/*" --ignore "*/all_multiqc/*"` so it no longer re-ingests MultiQC report/data dirs
- `nextflow.config` compacted with regex `withName` selectors (`fastqc.*`, `bbtools.*`, `ivar.*`,
  `samtools_(bam|coverage|mpileup)`); `samtools_screen` kept as its own block
- `modules/kraken2.nf` — reverted to bundled viral container (`staphb/kraken2:2.17.1-viral-20251015`); removed `--confidence` filter
- `nextflow.config` — reduced Kraken2 memory from 50 GB to 8 GB
- `modules/nextclade.nf` — switched to per-serotype community datasets (`community/v-gen-lab/dengue/denv1–4`)
- `bin/summary_report.py` — renamed `VADR_flag`/`QC_flag` → `vadr_flag`/`qc_flag`; removed `nextclade_qc_overall` column; improved unclassified/low-coverage QC message; VADR flag derived from QC flag for samples that never reached VADR
- `bin/qc_gate.py` — updated QC coverage threshold; depth fail message now shows actual observed depth
- `README.md` — corrected resource requirements; added conda setup step; added Nextflow version constraint note (≥23.04, <26.0); removed estimated runtime
- `daytona_dengue.sh` — updated SLURM resource allocation

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
