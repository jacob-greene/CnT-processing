/*
 * Genome preparation: canonical chromosome sizes, the effective genome size,
 * Bowtie2 indexes built from FASTA, and a GTF with UCSC chromosome names.
 * genome_size is always computed here from the chrom.sizes file, never hard-coded.
 */

process GENOME_SIZE {
    tag "${params.genome}"
    container params.container__bedtools
    publishDir "${params.outdir}/genome", mode: 'copy', pattern: '*.{sizes,tsv,bed}'

    input:
    path chrom_sizes
    path blacklist      // assets/NO_FILE when blacklist filtering is off
    val  mito_contig

    output:
    path 'canonical.chrom.sizes',   emit: chrom_sizes
    path 'genome_size.tsv',         emit: tsv
    path 'blacklist.canonical.bed', optional: true, emit: blacklist
    path 'versions.yml',            emit: versions

    script:
    def bl = blacklist.name != 'NO_FILE' ? blacklist : ''
    """
    genome_size.sh ${chrom_sizes} ${mito_contig} ${params.keep_chrM} ${bl}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}

/*
 * Build a Bowtie2 index from one or more FASTA files (plain or gzipped).
 * Used for E. coli (always, unless `ecoli_index` is given) and for the target
 * genome when `bowtie2_index = false`.
 * Output: a folder `bowtie2_<name>/` holding `<name>.*.bt2`.
 */
process BOWTIE2_BUILD {
    tag "${name}"
    label 'process_medium'
    container params.container__bowtie2

    input:
    val  name
    path fasta, stageAs: 'fasta/*'

    output:
    tuple path("bowtie2_${name}"), val(name), emit: index
    path 'versions.yml',                      emit: versions

    script:
    """
    mkdir bowtie2_${name}
    zcat -f fasta/* > ${name}.fa
    bowtie2-build --threads ${task.cpus} -q ${name}.fa bowtie2_${name}/${name}
    rm ${name}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | awk 'NR == 1 { print \$NF }')
    END_VERSIONS
    """
}

/*
 * GTF with UCSC chromosome names, gzipped.
 * With a chromosome map (two columns: GTF name, UCSC name), names are
 * translated and features on unlisted contigs are dropped. Without one, the
 * GTF is copied unchanged. Header lines (#...) are kept.
 */
process PREPARE_GTF {
    tag "${params.genome}"
    container params.container__bedtools
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path gtf
    path chrom_map      // assets/NO_FILE when the GTF already uses UCSC names

    output:
    path 'genes.gtf.gz', emit: gtf

    script:
    if (chrom_map.name == 'NO_FILE')
        """
        zcat -f ${gtf} | gzip -c > genes.gtf.gz
        """
    else
        """
        zcat -f ${gtf} \\
            | awk -F '\\t' -v OFS='\\t' '
                NR == FNR { if (\$0 !~ /^#/ && NF >= 2) map[\$1] = \$2; next }
                /^#/      { print; next }
                (\$1 in map) { \$1 = map[\$1]; print }
              ' ${chrom_map} - \\
            | gzip -c > genes.gtf.gz
        """
}
