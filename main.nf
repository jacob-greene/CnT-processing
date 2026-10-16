#!/usr/bin/env nextflow
/*
 * CnT-processing
 * CUT&Tag paired-end FASTQ -> coverage-normalized bigWigs, featureCounts matrices,
 * SEACR peaks, QC tables and a MultiQC report.
 *
 * Status: Phase 2. Implemented: sample sheet validation, genome preparation and
 * processing steps 1-10 (FastQC, trimming, alignment, E. coli QC, filtering,
 * duplicate marking, alignment stats, fragment sizes, fragments, bigWigs).
 * Steps 11-16 follow in Phase 3.
 */

nextflow.enable.dsl = 2

include { parseSamplesheet; samplesheetSummary } from './modules/samplesheet.nf'
include { GENOME_SIZE; PREPARE_GTF              } from './modules/prepare_genome.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_GENOME  } from './modules/prepare_genome.nf'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_ECOLI   } from './modules/prepare_genome.nf'
include { CAT_FASTQ; FASTQC; TRIM_GALORE         } from './modules/reads.nf'
include { BOWTIE2_ALIGN; BOWTIE2_ECOLI; FILTER_BAM; MARKDUPLICATES; NAMESORT } from './modules/align.nf'
include { FRAGMENTS; SCALE_FACTORS               } from './modules/coverage.nf'
include { FRAGMENTS as FRAGMENTS_DEDUP           } from './modules/coverage.nf'
include { BIGWIG as BIGWIG_COV                   } from './modules/coverage.nf'
include { BIGWIG as BIGWIG_RC                    } from './modules/coverage.nf'
include { QC_TABLES; FRAGMENT_SIZE_PLOT; SOFTWARE_VERSIONS } from './modules/qc.nf'

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
        fasta          : ref('fasta'),
    ]
    if (params.blacklist_filter && !genome.blacklist) {
        error "blacklist_filter = true, but no blacklist is configured for ${params.genome}. Set --blacklist."
    }

    def samples = parseSamplesheet(params.samplesheet)
    ch_versions = Channel.empty()

    log.info """\
        C n T - p r o c e s s i n g   ${workflow.manifest.version}
        ==========================================================
        samplesheet       : ${params.samplesheet}
        genome            : ${params.genome}
        outdir            : ${params.outdir}
        bowtie2_index     : ${genome.bowtie2_index ?: 'built from ' + genome.fasta}
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
    def no_file      = file("${projectDir}/assets/NO_FILE")
    def blacklist_in = params.blacklist_filter ? file(genome.blacklist) : no_file
    GENOME_SIZE(file(genome.chrom_sizes), blacklist_in, genome.mito_contig)
    ch_chrom_sizes = GENOME_SIZE.out.chrom_sizes
    ch_blacklist   = params.blacklist_filter ? GENOME_SIZE.out.blacklist : Channel.value(no_file)

    GENOME_SIZE.out.tsv
        .map { tsv -> tsv.readLines().collectEntries { it.split('\t') as List } }
        .subscribe { gs -> log.info "genome_size used: ${gs.genome_size} bp (blacklist_filter = ${gs.blacklist_filter})" }

    // GTF with UCSC chromosome names (used by featureCounts and the TSS plots).
    PREPARE_GTF(file(genome.gtf), genome.gtf_chrom_map ? file(genome.gtf_chrom_map) : no_file)

    // Bowtie2 indexes: the target genome (given, or built from `fasta`) and E. coli.
    if (genome.bowtie2_index) {
        def prefix = file(genome.bowtie2_index)
        ch_index = Channel.value([prefix.parent, prefix.name])
    } else {
        if (!genome.fasta) {
            error "bowtie2_index = false, but no --fasta is given to build the index from."
        }
        def fastas = (genome.fasta instanceof List ? genome.fasta : [genome.fasta]).collect { file(it) }
        BOWTIE2_BUILD_GENOME(params.genome, fastas)
        ch_index = BOWTIE2_BUILD_GENOME.out.index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD_GENOME.out.versions)
    }
    if (params.ecoli_index) {
        def prefix = file(params.ecoli_index)
        ch_ecoli_index = Channel.value([prefix.parent, prefix.name])
    } else {
        BOWTIE2_BUILD_ECOLI('ecoli', [file(params.ecoli_fasta)])
        ch_ecoli_index = BOWTIE2_BUILD_ECOLI.out.index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD_ECOLI.out.versions)
    }

    // ---- Steps 1-2: lanes, FastQC, trimming ----
    ch_reads = Channel.fromList(samples).map { s -> [s.meta, s.fastq_1, s.fastq_2] }
    CAT_FASTQ(ch_reads)
    FASTQC(CAT_FASTQ.out.reads)
    TRIM_GALORE(CAT_FASTQ.out.reads)

    // ---- Steps 3-4: alignment to the target genome and to E. coli ----
    BOWTIE2_ALIGN(TRIM_GALORE.out.reads, ch_index)
    BOWTIE2_ECOLI(TRIM_GALORE.out.reads, ch_ecoli_index)

    // ---- Step 5: pair filter, sort, index; % chrM and fragment lengths ----
    FILTER_BAM(BOWTIE2_ALIGN.out.bam, genome.mito_contig)

    // ---- Step 6: mark duplicates (kept) ----
    ch_markdup_in = FILTER_BAM.out.bam
        .join(FILTER_BAM.out.stats)
        .map { meta, bam, bai, stats ->
            def illumina = stats.readLines().find { it.startsWith('illumina_read_names\t') }.endsWith('true')
            [meta, bam, bai, illumina]
        }
    MARKDUPLICATES(ch_markdup_in)
    MARKDUPLICATES.out.metrics
        .filter { meta, metrics, illumina -> !illumina }
        .collect { meta, metrics, illumina -> meta.id }
        .subscribe { ids ->
            log.warn "Read names have no Illumina flow-cell coordinates (e.g. SRA-renamed reads) for: " +
                     "${ids.join(', ')}. Optical duplicates are reported as NA."
        }

    // ---- Step 9: fragments (+ optional blacklist filter) and raw bedGraphs ----
    NAMESORT(MARKDUPLICATES.out.bam)
    FRAGMENTS(NAMESORT.out.bam, ch_chrom_sizes, ch_blacklist, '')
    // dedup_controls = true: deduplicated control fragments + bedGraph, used only by SEACR.
    FRAGMENTS_DEDUP(NAMESORT.out.dedup, ch_chrom_sizes, ch_blacklist, '.dedup')

    // ---- Step 10: scale factors and bigWigs ----
    SCALE_FACTORS(
        Channel.value(samples.collect { "${it.meta.id}=${it.meta.pool}" }),
        FRAGMENTS.out.stats.map { meta, f -> f }.collect(),
        GENOME_SIZE.out.tsv
    )
    ch_factors = SCALE_FACTORS.out.tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> [row.sample, row.s_cov, row.s_rc] }
    ch_bg_factors = FRAGMENTS.out.bedgraph
        .map { meta, bg -> [meta.id, meta, bg] }
        .join(ch_factors)
        .map { id, meta, bg, s_cov, s_rc -> [meta, bg, s_cov, s_rc] }

    BIGWIG_COV(ch_bg_factors.map { meta, bg, s_cov, s_rc -> [meta, bg, s_cov] }, ch_chrom_sizes, 'coverage_norm')
    if (params.readcount_scaled_bw) {
        BIGWIG_RC(ch_bg_factors.map { meta, bg, s_cov, s_rc -> [meta, bg, s_rc] }, ch_chrom_sizes, 'readcount_scaled')
    }

    // ---- Steps 7-8: summary tables and the fragment size plot ----
    ch_qc_files = CAT_FASTQ.out.stats
        .mix(BOWTIE2_ALIGN.out.log, BOWTIE2_ECOLI.out.log, FILTER_BAM.out.stats,
             MARKDUPLICATES.out.metrics.map { meta, m, illumina -> [meta, m] }, FRAGMENTS.out.stats)
        .map { meta, f -> f }
        .collect()
    QC_TABLES(Channel.value(samples.collect { it.meta.id }), ch_qc_files)
    FRAGMENT_SIZE_PLOT(FILTER_BAM.out.fragment_sizes.map { meta, f -> f }.collect())

    // ---- software_versions.yml ----
    ch_versions = ch_versions.mix(
        GENOME_SIZE.out.versions, FASTQC.out.versions, TRIM_GALORE.out.versions,
        BOWTIE2_ALIGN.out.versions, FILTER_BAM.out.versions, MARKDUPLICATES.out.versions,
        FRAGMENTS.out.versions, BIGWIG_COV.out.versions)
    SOFTWARE_VERSIONS(ch_versions.collect())

    // Steps 11-16 (featureCounts, SEACR, FRiP, deepTools, correlation, MultiQC) follow in Phase 3.
}
