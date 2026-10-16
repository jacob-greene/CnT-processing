#!/usr/bin/env nextflow
/*
 * CnT-processing
 * CUT&Tag paired-end FASTQ -> coverage-normalized bigWigs, featureCounts matrices,
 * SEACR peaks, QC tables and a MultiQC report.
 *
 * Status: Phase 1 scaffold. Implemented: sample sheet validation, genome
 * resolution and genome_size computation. Processing steps 1-16 follow.
 */

nextflow.enable.dsl = 2

include { parseSamplesheet; samplesheetSummary } from './modules/samplesheet.nf'
include { GENOME_SIZE                          } from './modules/prepare_genome.nf'

def SUPPORTED_GENOMES = ['hg38', 'hg19', 'mm10', 'dm6']

// A reference value: the command-line / params.json override, else the genomes.config entry.
// An empty string (e.g. a blank BASH Workbench field) counts as "not set".
def ref(String key) {
    def override = params[key]
    return (override != null && override != '') ? override : params.genomes[params.genome][key]
}

workflow {

    // ---- Parameter checks: fail before any job starts ----
    if (!params.samplesheet) {
        error "Missing required parameter: --samplesheet <file.csv>"
    }
    if (!params.genome) {
        error "Missing required parameter: --genome (one of ${SUPPORTED_GENOMES.join(', ')})"
    }
    if (!(params.genome in SUPPORTED_GENOMES)) {
        error "Unknown genome '${params.genome}'. Supported genomes: ${SUPPORTED_GENOMES.join(', ')}"
    }

    def genome = [
        chrom_sizes    : ref('chrom_sizes'),
        mito_contig    : ref('mito_contig'),
        gtf            : ref('gtf'),
        gtf_chrom_map  : ref('gtf_chrom_map'),
        ccre_bed       : ref('ccre_bed'),
        blacklist      : ref('blacklist'),
        bowtie2_index  : ref('bowtie2_index'),
    ]
    if (params.blacklist_filter && !genome.blacklist) {
        error "blacklist_filter = true, but no blacklist is configured for ${params.genome}. Set --blacklist."
    }

    def samples = parseSamplesheet(params.samplesheet)

    log.info """\
        C n T - p r o c e s s i n g   ${workflow.manifest.version}
        ==========================================================
        samplesheet       : ${params.samplesheet}
        genome            : ${params.genome}
        outdir            : ${params.outdir}
        bowtie2_index     : ${genome.bowtie2_index}
        chrom_sizes       : ${genome.chrom_sizes}
        gtf               : ${genome.gtf}
        ccre_bed          : ${genome.ccre_bed ?: 'none (cCRE counting skipped)'}
        blacklist_filter  : ${params.blacklist_filter}${params.blacklist_filter ? ' (' + genome.blacklist + ')' : ''}
        keep_chrM         : ${params.keep_chrM}
        min_mapq          : ${params.min_mapq}
        readcount_scaled  : ${params.readcount_scaled_bw}
        dedup_controls    : ${params.dedup_controls}

        Samples (${samples.size()}):
        ${samplesheetSummary(samples).replace('\n', '\n        ')}
        """.stripIndent()

    if (!genome.ccre_bed) {
        log.info "No ENCODE cCREs for ${params.genome}: cCRE counting will be skipped."
    }

    if (params.validate_only) {
        log.info "validate_only = true: the sample sheet and genome settings are valid. Stopping."
        return
    }

    // ---- Genome preparation ----
    def blacklist_in = params.blacklist_filter ? file(genome.blacklist) : file("${projectDir}/assets/NO_FILE")
    GENOME_SIZE(file(genome.chrom_sizes), blacklist_in, genome.mito_contig)

    GENOME_SIZE.out.tsv
        .map { tsv -> tsv.readLines().collectEntries { it.split('\t') as List } }
        .subscribe { gs -> log.info "genome_size used: ${gs.genome_size} bp (blacklist_filter = ${gs.blacklist_filter})" }

    // ---- Reads ----
    ch_reads = Channel.fromList(samples).map { s -> [s.meta, s.fastq_1, s.fastq_2] }

    // Processing steps 1-16 are added in later phases.
}
