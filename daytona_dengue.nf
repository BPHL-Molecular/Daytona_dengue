#!/usr/bin/env nextflow

/*
  Daytona Dengue
  Florida's BPHL Nextflow pipeline for Dengue virus NGS data analysis
  Authors: Yibo Dong, Arnold Rodriguez
  Email: bphl-sebioinformatics@flhealth.gov
*/

nextflow.enable.dsl = 2

params.assets_dir = "${projectDir}/assets"

include { fastqc; fastqc_clean             } from './modules/fastqc.nf'
include { humanscrubber                    } from './modules/humanscrubber.nf'
include { trimmomatic                      } from './modules/trimmomatic.nf'
include { bbtools_adapters; bbtools_phix   } from './modules/bbtools.nf'
include { multiqc                          } from './modules/multiqc.nf'
include { kraken2                          } from './modules/kraken2.nf'
include { bwa                              } from './modules/bwa.nf'
include { serotype_detect                  } from './modules/serotype_detect.nf'
include { samtools_screen                  } from './modules/samtools.nf'
include { samtools_bam                     } from './modules/samtools.nf'
include { samtools_coverage                } from './modules/samtools.nf'
include { samtools_mpileup                 } from './modules/samtools.nf'
include { ivar_trim                        } from './modules/ivar.nf'
include { ivar_variants                    } from './modules/ivar.nf'
include { ivar_consensus                   } from './modules/ivar.nf'
include { qc_gate                          } from './modules/qc_gate.nf'
include { vadr_download; vadr              } from './modules/vadr.nf'
include { nextclade_download; nextclade    } from './modules/nextclade.nf'
include { summary_report                   } from './modules/summary_report.nf'

def refFileList(String sero) {
    def nameMap = [
        DENV1: 'NC_001477.1_DENV1',
        DENV2: 'NC_001474.2_DENV2',
        DENV3: 'NC_001475.2_DENV3',
        DENV4: 'NC_002640.1_DENV4',
    ]
    def fa = "${projectDir}/assets/reference/${nameMap[sero]}.fasta"
    return [fa, "${fa}.amb", "${fa}.ann", "${fa}.bwt",
            "${fa}.fai", "${fa}.pac", "${fa}.sa"].collect { p -> file(p) }
}

workflow {

    log.info """
    Daytona Dengue — BPHL Dengue NGS pipeline
    ==========================================================================
    input dir    : ${params.input}
    output dir   : ${params.output}
    assets       : ${projectDir}/assets
    ==========================================================================
    """

    ch_reads = channel.fromFilePairs(
        ["${params.input}/*_{1,2}.fastq.gz",
         "${params.input}/*_R{1,2}_*.fastq.gz"],
        checkIfExists: false
    )
    .map { id, files ->
        def clean_id = id.replaceAll(/_S\d+_L\d+$/, '')
        def meta     = [ id: clean_id, single_end: false ]
        [ meta, files ]
    }

    ch_sd_denv1 = channel.value(refFileList('DENV1'))
    ch_sd_denv2 = channel.value(refFileList('DENV2'))
    ch_sd_denv3 = channel.value(refFileList('DENV3'))
    ch_sd_denv4 = channel.value(refFileList('DENV4'))

    fastqc(ch_reads)
    humanscrubber(ch_reads)
    trimmomatic(humanscrubber.out.reads)
    bbtools_adapters(trimmomatic.out.reads)
    bbtools_phix(bbtools_adapters.out.reads)
    fastqc_clean(bbtools_phix.out.reads)

    ch_multiqc_input = fastqc.out.zip
        .mix(fastqc_clean.out.zip)
        .map { _meta, zip -> zip }
        .collect()
    multiqc(ch_multiqc_input)

    kraken2(bbtools_phix.out.reads)

    bwa(bbtools_phix.out.reads, ch_sd_denv1, ch_sd_denv2, ch_sd_denv3, ch_sd_denv4)
    samtools_screen(bwa.out.sams)
    serotype_detect(samtools_screen.out.coverage)

    ch_serotyped = serotype_detect.out.serotype
        .map { meta, serotype_file ->
            def serotype = serotype_file.text.trim()
            [ meta + [ serotype: serotype ], meta ]
        }
        .branch { enriched_meta, _orig ->
            typed:   enriched_meta.serotype != 'unclassified'
            untyped: enriched_meta.serotype == 'unclassified'
        }

    ch_serotyped.untyped
        .map { enriched_meta, _orig -> enriched_meta.id }
        .collect()
        .subscribe { ids ->
            if (ids) log.warn "Unserotyped samples (excluded from pipeline): ${ids.join(', ')}"
        }

    ch_clean_typed = ch_serotyped.typed
        .map { enriched_meta, orig_meta -> [ orig_meta.id, enriched_meta ] }
        .join(
            bbtools_phix.out.reads.map { meta, reads -> [ meta.id, reads ] }
        )
        .map { _id, enriched_meta, reads -> [ enriched_meta, reads ] }

    def primer_map = [
        DENV1: "${projectDir}/assets/primers/DENV1.primer.bed",
        DENV2: "${projectDir}/assets/primers/DENV2.primer.bed",
        DENV3: "${projectDir}/assets/primers/DENV3.primer.bed",
        DENV4: "${projectDir}/assets/primers/DENV4.primer.bed",
    ]
    def gff_map = [
        DENV1: "${projectDir}/assets/annotations/Dengue1_GCF_000862125.1_ViralProj15306_genomic.gff.gz",
        DENV2: "${projectDir}/assets/annotations/Dengue2_GCF_000871845.1_ViralProj20183_genomic.gff.gz",
        DENV3: "${projectDir}/assets/annotations/Dengue3_GCF_000866625.1_ViralProj15598_genomic.gff.gz",
        DENV4: "${projectDir}/assets/annotations/Dengue4_GCF_000865065.1_ViralProj15599_genomic.gff.gz",
    ]

    ch_clean_with_ref = ch_clean_typed
        .map { meta, _reads ->
            [ meta, refFileList(meta.serotype),
              file(primer_map[meta.serotype]),
              file(gff_map[meta.serotype]) ]
        }

    ch_winning_sam = ch_clean_typed
        .map { meta, _reads -> [ meta.id, meta ] }
        .join(
            bwa.out.sams.map { meta, sam_files ->
                def files = sam_files instanceof List ? sam_files : [sam_files]
                [ meta.id, files ]
            }
        )
        .map { _id, meta, sam_files ->
            def winning = sam_files.find { f -> f.name.contains("_${meta.serotype}.sam") }
            [ meta, winning ]
        }

    samtools_bam(ch_winning_sam)

    ch_bam_with_primer = samtools_bam.out.bam
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join(
            ch_clean_with_ref.map { meta, _ref_list, primer, _gff -> [ meta.id, primer ] }
        )
        .map { _id, meta, bam, bai, primer -> [ meta, bam, bai, primer ] }

    ivar_trim(ch_bam_with_primer)

    samtools_coverage(ivar_trim.out.bam)

    ch_trimmed_bam_with_ref = ivar_trim.out.bam
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join(
            ch_clean_with_ref.map { meta, ref_list, _primer, _gff ->
                def ref_fa = ref_list instanceof List ? ref_list.find { f -> f.name.endsWith('.fasta') } : ref_list
                [ meta.id, ref_fa ]
            }
        )
        .map { _id, meta, bam, bai, ref -> [ meta, bam, bai, ref ] }

    samtools_mpileup(ch_trimmed_bam_with_ref)

    ch_mpileup_with_ref_gff = samtools_mpileup.out.mpileup
        .map { meta, mpileup -> [ meta.id, meta, mpileup ] }
        .join(
            ch_clean_with_ref.map { meta, ref_list, _primer, gff ->
                def ref_fa = ref_list instanceof List ? ref_list.find { f -> f.name.endsWith('.fasta') } : ref_list
                [ meta.id, ref_fa, gff ]
            }
        )
        .map { _id, meta, mpileup, ref, gff -> [ meta, mpileup, ref, gff ] }

    ivar_variants(ch_mpileup_with_ref_gff)

    ivar_consensus(samtools_mpileup.out.mpileup)

    ch_qc_input = ivar_consensus.out.consensus
        .join(samtools_coverage.out.coverage)

    qc_gate(ch_qc_input)

    ch_qc_pass = qc_gate.out.qc
        .filter { _meta, qc_file ->
            def lines = qc_file.readLines()
            lines.size() > 1 && lines[1].split('\t')[1].trim() == 'PASS'
        }

    ch_vadr_input = ch_qc_pass
        .map { meta, _qc -> [ meta.id, meta ] }
        .join( ivar_consensus.out.consensus.map { meta, f -> [ meta.id, f ] } )
        .map { _id, meta, consensus -> [ meta, consensus ] }

    ch_vadr_models = vadr_download()
    vadr(ch_vadr_input, ch_vadr_models.models)

    nextclade_download(channel.of('DENV1', 'DENV2', 'DENV3', 'DENV4'))

    ch_nextclade_input = ivar_consensus.out.consensus
        .map { meta, consensus -> [meta.serotype, meta, consensus] }
        .combine(nextclade_download.out.db, by: 0)
        .map { _sero, meta, consensus, dataset -> [meta, consensus, dataset] }

    nextclade(ch_nextclade_input)

    ch_barrier = vadr.out.done
        .mix(nextclade.out.done)
        .map { _meta -> 1 }
        .collect()
        .map { _ids -> true }

    summary_report(
        ch_barrier,
        qc_gate.out.qc.map                                       { _meta, f -> f }.collect().ifEmpty([]),
        samtools_coverage.out.coverage.map                       { _meta, f -> f }.collect().ifEmpty([]),
        ivar_consensus.out.consensus.map                         { _meta, f -> f }.collect().ifEmpty([]),
        nextclade.out.tsv.map                                    { _meta, f -> f }.collect().ifEmpty([]),
        vadr.out.results.map                                     { _meta, f -> f }.collect().ifEmpty([]),
        kraken2.out.report.map                                   { _meta, f -> f }.collect(),
        serotype_detect.out.serotype.map                         { _meta, f -> f }.collect(),
        samtools_screen.out.coverage.flatMap                     { _meta, fs -> fs }.collect(),
        trimmomatic.out.stats.map                                { _meta, f -> f }.collect(),
        bbtools_phix.out.phix_log.map                            { _meta, f -> f }.collect()
    )
}
