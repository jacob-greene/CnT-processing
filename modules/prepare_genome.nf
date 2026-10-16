/*
 * Genome preparation: canonical chromosome sizes and the effective genome size.
 * genome_size is always computed here from the chrom.sizes file, never hard-coded.
 */

process GENOME_SIZE {
    tag "${params.genome}"
    container params.container__bedtools
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path chrom_sizes
    path blacklist      // assets/NO_FILE when blacklist filtering is off
    val  mito_contig

    output:
    path 'canonical.chrom.sizes', emit: chrom_sizes
    path 'genome_size.tsv',       emit: tsv
    path 'blacklist.canonical.bed', optional: true, emit: blacklist

    script:
    def bl = blacklist.name != 'NO_FILE' ? blacklist : ''
    """
    genome_size.sh ${chrom_sizes} ${mito_contig} ${params.keep_chrM} ${bl}
    """
}
