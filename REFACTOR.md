# Daytona Dengue — Refactor Tracking

BPHL Sanibel V2.1 modernization + Nextclade integration.  
**Key decisions:** `meta.serotype` drives all serotype routing (no more `if/elif SER1/2/3/4`); Kraken2 runs inside the pipeline; Nextclade (`nextstrain/dengue/all`) adds clade assignment; VADR kept for GenBank submission; all post-processing scripts retired.

---

## Phase 1 — Repository Restructuring

- [x] Create `bin/` directory
- [x] Create `assets/annotations/` directory
- [x] Move GFF files from `reference/` → `assets/annotations/`
  - `Dengue1_GCF_000862125.1_ViralProj15306_genomic.gff.gz`
  - `Dengue2_GCF_000871845.1_ViralProj20183_genomic.gff.gz`
  - `Dengue3_GCF_000866625.1_ViralProj15598_genomic.gff.gz`
  - `Dengue4_GCF_000865065.1_ViralProj15599_genomic.gff.gz`
- [x] Delete `configs/` directory (completed end of Phase 2)
- [x] Add `CHANGELOG.md` to repo root
- [x] Update `README.md` to reflect new structure

**Retired files (deleted):**
- `kraken2_viral.sh` — replaced by `modules/kraken2.nf` + `bin/parse_kraken2.py`
- `renamefile.sh` — replaced by `id.replaceAll()` in input channel
- `table.py` — replaced by `bin/summary_report.py`
- `rename_aqp.py` — replaced by meta-driven design (no SER prefix in filenames)

---

## Phase 2 — `nextflow.config` (full rewrite)

> Replace the current 3-line stub and `configs/` with a single complete config.

- [x] No `params {}` block — tool parameters live in their respective modules; pipeline paths live only in `params.yaml`
- [x] `withName:` block for every process (container + cpus + memory):
  - `kraken2` — `docker://staphb/kraken2:2.1.3`, cpus=8, mem=32.GB
  - `fastqc` — `docker://staphb/fastqc:0.11.9`, cpus=2, mem=4.GB
  - `fastqc_clean` — `docker://staphb/fastqc:0.11.9`, cpus=2, mem=4.GB
  - `humanscrubber` — `docker://ncbi/sra-human-scrubber:1.1.2021-05-05`, cpus=4, mem=8.GB
  - `trimmomatic` — `docker://staphb/trimmomatic:0.39`, cpus=4, mem=8.GB
  - `bbtools_adapters` — `docker://staphb/bbtools:38.76`, cpus=4, mem=8.GB
  - `bbtools_phix` — `docker://staphb/bbtools:38.76`, cpus=4, mem=8.GB
  - `multiqc` — `docker://staphb/multiqc:1.8`, cpus=2, mem=4.GB
  - `bwa_mem` — `docker://staphb/bwa:0.7.17`, cpus=8, mem=16.GB
  - `samtools_bam`, `samtools_coverage`, `samtools_mpileup` — `docker://staphb/samtools:1.12`, cpus=4, mem=8.GB
  - `ivar_trim`, `ivar_variants`, `ivar_consensus` — `docker://staphb/ivar:1.3.2`, cpus=4, mem=8.GB
  - `pystats` — no container (bin/ script via conda), cpus=1, mem=4.GB
  - `vadr` — `docker://staphb/vadr:1.5.1`, cpus=8, mem=16.GB, `errorStrategy = 'ignore'`
  - `nextclade_download` — `docker://nextstrain/nextclade:3.0.0`, cpus=2, mem=4.GB
  - `nextclade` — `docker://nextstrain/nextclade:3.0.0`, cpus=4, mem=8.GB
  - `summary_report` — no container (bin/ script via conda), cpus=1, mem=2.GB
- [x] Profiles: `standard`, `docker`, `singularity`, `apptainer`
- [x] Global `process.errorStrategy = 'finish'` and `process.maxRetries = 0` at bottom
- [x] No `manifest` block, no `check_max()`, no `includeConfig` calls

---

## Phase 3 — `params.yaml` (minimal update)

> Only runtime paths that have no sensible default live here. All tool parameters
> (Nextclade dataset name/tag/cache, Kraken2 confidence, trimming thresholds, iVar
> settings, QC cutoffs) are defaults inside their respective modules.

- [x] `input:` and `output:` listed first, no trailing slash
- [x] Add `reference:`, `primer:`, `annotations:` keys
- [x] Add `kraken_db:` (required — large external DB, no bundled default)
- [x] Remove HPC-specific absolute paths from original file

---

## Phase 4 — `daytona_dengue.sh` (SLURM script update)

- [x] Update SBATCH header: `--job-name=daytona_dengue`, add `--mail-user` and `--mail-type=FAIL,END`
- [x] Load modules: `module load conda nextflow apptainer`
- [x] Activate environment: `conda activate PIPELINE_ENV`
- [x] Export `NXF_APPTAINER_CACHEDIR`
- [x] Run command: `nextflow run daytona_dengue.nf -profile apptainer -params-file params.yaml`
- [x] Add timestamp rename block (on success: `mv output output-YYYYMMDDHHMMSS`)
- [x] Add commented-out `rm -rf ./work ./cache` cleanup line
- [x] **Remove** all post-processing: `sort`, `mkdir`, `mv`, `python3 table.py`, `python3 rename_aqp.py`, `rm -r ./output/dengue*`

---

## Phase 5 — `daytona_dengue.nf` (full rewrite)

- [x] DSL2 header + comment block (pipeline name, purpose, authors, email)
- [x] `nextflow.enable.dsl = 2`
- [x] `include {}` for all modules (one per line, lowercase process names)
- [x] `log.info` block printing `input`, `output`, `kraken_db`, `nextclade_dataset_name`
- [x] Input channel via `channel.fromFilePairs` — both `*_{1,2}.fastq.gz` and `*_R{1,2}_*.fastq.gz`; `checkIfExists: false`
- [x] `id.replaceAll(/_S\d+_L\d+$/, '')` to strip Illumina suffix
- [x] `meta = [ id: clean_id, single_end: false ]`
- [x] `kraken2(ch_reads)` → `bin/parse_kraken2.py` → enrich meta with `meta.serotype`
- [x] Unserotyped branch: `.branch {}` → `ch_unserotyped` logs warning; `ch_typed` continues
- [x] Serotype lookup maps (inline `.map {}`) for `reference`, `primer`, `gff` paths
- [x] All channel variables prefixed with `ch_`
- [x] Pipeline chain (named `ch_` variables, no pipe operator):
  `humanscrubber → trimmomatic → bbtools_adapters → bbtools_phix → fastqc + fastqc_clean → multiqc → bwa_mem → samtools_bam → ivar_trim → samtools_coverage + samtools_mpileup → ivar_variants + ivar_consensus → pystats`
- [x] QC filter: `ch_qc_pass = pystats.out.filter { meta, _s -> meta.qc_pass }`
- [x] `nextclade_download()` → `ch_nextclade_db` (storeDir, runs once)
- [x] `nextclade(ivar_consensus.out.consensus, ch_nextclade_db)` — all consensus sequences
- [x] `vadr(ch_qc_pass)` — QC-pass samples only
- [x] Barrier channel: `mix(vadr.out.done, nextclade.out.done).collect()` before `summary_report`
- [x] `summary_report(barrier, pystats_files, nextclade_tsvs, vadr_results)`

---

## Phase 6 — Module Rewrites

### Group A (independent — can be written in parallel)

#### `modules/kraken2.nf` ⬅ NEW
- [x] Process `kraken2`: input `tuple val(meta), path(reads)`; outputs report + serotype files
- [x] Runs `kraken2 --db ${params.kraken_db} --confidence 0.5 --paired`
- [x] Calls `bin/parse_kraken2.py` to extract serotype from report
- [x] `publishDir "${params.output}/${meta.id}/kraken2", mode: 'copy'`
- [x] `emit: report`, `emit: serotype`

#### `modules/fastqc.nf` ⬅ modernize + absorb `fastqc_clean.nf`
- [x] Process `fastqc`: raw reads QC; `tag`, `publishDir`, `emit: zip`
- [x] Process `fastqc_clean`: clean reads QC (separate process, same file)
- [x] Remove `mkdir` logic; use `task.cpus`
- [x] Delete old `modules/fastqc_clean.nf`

#### `modules/humanscrubber.nf` ⬅ modernize
- [x] Add `tag "${meta.id}"`, `publishDir "${params.output}/${meta.id}/humanscrubber", mode: 'copy'`
- [x] Tuple input/output: `tuple val(meta), path(reads)`
- [x] Named `emit: reads`; use `task.cpus`; remove hardcoded paths

#### `modules/trimmomatic.nf` ⬅ modernize
- [x] Hardcoded defaults in module: `SLIDINGWINDOW:4:30 MINLEN:75 TRAILING:20`
- [x] Use `task.cpus`; remove if/elif SER1–4 blocks; named `emit: reads`

---

### Group B (independent — can be written in parallel)

#### `modules/bbtools.nf` ⬅ rename from `bbduk.nf`
- [x] Process `bbtools_adapters`: adapter trimming, `emit: reads`
- [x] Process `bbtools_phix`: PhiX removal, `emit: reads`
- [x] Both: `tag`, `publishDir`, `task.cpus`; remove if/elif SER1–4 blocks
- [x] Delete old `modules/bbduk.nf`

#### `modules/multiqc.nf` ⬅ modernize
- [x] Input: collected fastqc zips + trimmomatic logs
- [x] `publishDir "${params.output}/multiqc", mode: 'copy'`; `emit: report`
- [x] Remove placeholder directory creation

#### `modules/bwa.nf` ⬅ rename from `frag_bwa.nf`
- [x] Process `bwa_mem`: input `tuple val(meta), path(reads), path(reference)`
- [x] Reference as channel input — no if/elif; `emit: sam`; `task.cpus`
- [x] Delete old `modules/frag_bwa.nf`

---

### Group C (depends on Group A + B)

#### `modules/samtools.nf` ⬅ merge `frag_samtools.nf` + `primer_trim_samtools.nf`
- [x] Process `samtools_bam`: SAM → sorted + indexed BAM; `emit: bam`
- [x] Process `samtools_coverage`: coverage stats from trimmed BAM; `emit: coverage`
- [x] Process `samtools_mpileup`: pileup; uses hardcoded `-d 8000 -Q 0`; reference passed as channel input; `emit: mpileup`
- [x] All: `tag`, `publishDir`, `task.cpus`; remove if/elif blocks; remove redundant intermediate sorts
- [x] Delete old `modules/frag_samtools.nf` and `modules/primer_trim_samtools.nf`

#### `modules/ivar.nf` ⬅ merge `primer_trim_ivar.nf` + `assembly.nf`
- [x] Process `ivar_trim`: input `tuple val(meta), path(bam), path(bai), path(primer)`; sorts + indexes output; `emit: bam`
- [x] Process `ivar_variants`: input `tuple val(meta), path(mpileup), path(reference), path(gff)`; hardcoded `-m 10 -q 20 -t 0.25`; `emit: variants`
- [x] Process `ivar_consensus`: input `tuple val(meta), path(mpileup)`; hardcoded `-t 0 -m 10 -n N`; `emit: consensus`
- [x] All: `tag`, `publishDir`, `task.cpus`; no if/elif; no hardcoded GFF filenames
- [x] Delete old `modules/primer_trim_ivar.nf` and `modules/assembly.nf`

---

### Group D (depends on Group C)

#### `modules/pystats.nf` ⬅ thin wrapper only
- [x] Remove all embedded Python (~800 lines)
- [x] Script block calls `bin/pystats.py` with explicit CLI args
- [x] Input: `tuple val(meta), path(consensus), path(coverage), path(variants)`
- [x] Output: `tuple val(meta), path("${prefix}_stats.tsv"), emit: stats`
- [x] `publishDir "${params.output}/${meta.id}/stats", mode: 'copy'`

#### `modules/vadr.nf` ⬅ NEW
- [x] Process `vadr`: input `tuple val(meta), path(consensus)`
- [x] Runs `fasta-trim-terminal-ambigs.pl` then `v-annotate.pl --group Dengue --mkey flavi --mdir /opt/vadr/vadr-models-flavi/`
- [x] `errorStrategy = 'ignore'` declared in `nextflow.config` withName block
- [x] `emit: results` (path to vadr_results/) and `emit: done` (val meta, for barrier)
- [x] `publishDir "${params.output}/${meta.id}/vadr", mode: 'copy'`

#### `modules/nextclade_download.nf` ⬅ NEW
- [x] `storeDir params.output/db/nextclade` — skips re-download on resume
- [x] `nextclade dataset get --name 'nextstrain/dengue/all'`
- [x] `emit: db`

#### `modules/nextclade.nf` ⬅ NEW
- [x] Input: `tuple val(meta), path(consensus)` + `path(dataset)`
- [x] `nextclade run --input-dataset ${dataset} --output-tsv ${prefix}_nextclade.tsv ${consensus}`
- [x] `emit: tsv` and `emit: done` (val meta, for barrier)
- [x] `publishDir "${params.output}/${meta.id}/nextclade", mode: 'copy'`

#### `modules/summary_report.nf` ⬅ NEW
- [x] Input: `val barrier` + collected stats / nextclade / vadr paths
- [x] Calls `bin/summary_report.py`
- [x] `publishDir "${params.output}", mode: 'copy'`; `emit: report`

---

## Phase 7 — `bin/` Scripts

#### `bin/parse_kraken2.py` ⬅ NEW
- [x] `#!/usr/bin/env python3`; `chmod +x`
- [x] Reads Kraken2 report; identifies DENV1–4 from taxon name
- [x] Checks confidence >= threshold; outputs one-line `DENV1` (or `unclassified`)
- [x] Exits non-zero on parse failure

#### `bin/pystats.py` ⬅ NEW (extracted + cleaned from `pystats.nf`)
- [x] `#!/usr/bin/env python3`; `chmod +x`
- [x] CLI args: `--consensus`, `--coverage`, `--variants`, `--sample-id`, `--serotype`, `--qc-min-coverage`, `--qc-min-depth`
- [x] No Nextflow variable interpolation; no if/elif SER1/2/3/4; no subprocess for file counting
- [x] No VADR calls — VADR is now a separate Nextflow process
- [x] Outputs `{sample_id}_stats.tsv`

#### `bin/summary_report.py` ⬅ NEW
- [x] `#!/usr/bin/env python3`; `chmod +x`
- [x] CLI args: `--stats-dir`, `--nextclade-dir`, `--vadr-dir`, `--output`
- [x] Reads `*_stats.tsv`, `*_nextclade.tsv`, VADR `pass.list`/`fail.list`; joins on sample ID
- [x] Outputs `summary_report.tsv` columns: `sample_id`, `serotype`, `mapped_pct`, `genome_coverage_pct`, `mean_depth`, `n_bases`, `qc_flag`, `nextclade_clade`, `nextclade_qc`, `vadr_result`

---

## Verification Checklist

- [ ] `nextflow run daytona_dengue.nf -dry-run -params-file params.yaml` — no parse errors
- [ ] Channel graph trace confirms `meta.serotype` present on all downstream processes
- [ ] No `if`/`elif` blocks referencing serotype in any `.nf` module file
- [ ] No hardcoded absolute paths in any `.nf` file
- [ ] End-to-end test with ≥2 serotypes (e.g., DENV1 + DENV2 samples)
- [ ] Per-sample output directories created under `params.output/{sample_id}/`
- [ ] Nextclade dataset cached in `db/nextclade/` — second run skips download
- [ ] VADR runs only on QC-pass samples; QC-fail samples skipped without pipeline crash
- [ ] `summary_report.tsv` at `params.output/` root with all expected columns
- [ ] Unserotyped samples logged to warning file; pipeline does not crash
- [ ] `-resume` works correctly after a partial run

---

## Notes

- **Nextclade dataset:** `nextstrain/dengue/all` — covers DENV1–4 via phylogenetic placement; no serotype-specific dataset needed
- **GFF files:** Used by `ivar_variants -g` for codon-aware variant annotation. Nextclade uses its own bundled `genome_annotation.gff3` internally — separate files, separate tools
- **VADR model path:** `/opt/vadr/vadr-models-flavi/` is bundled inside `staphb/vadr:1.5.1` — no external path needed
- **`ivar` container:** Original config had no version tag on `primer_trim_ivar` — pinned to `1.3.2` in new config
- **Container to verify before implementation:** `docker://nextstrain/nextclade:3.0.0` — confirm latest available tag
