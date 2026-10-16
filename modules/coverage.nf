/*
 * Steps 9-10: fragment BED files, optional blacklist filtering, scale factors,
 * raw-count bedGraphs and normalized bigWigs.
 *
 * Normalization (see README, "Normalization"):
 *
 *   C_i         = total per-base coverage of sample i = sum of its fragment lengths
 *   genome_size = sum of canonical chromosome sizes (minus blacklisted bp when
 *                 blacklist_filter = true)
 *
 *   (a) coverage-normalized:  s_cov_i = genome_size / C_i
 *       Each base holds (fraction of total coverage at that base) x genome_size,
 *       so a uniform coverage would be 1 everywhere, and the genome-length-weighted
 *       mean of every bigWig is 1.
 *
 *   (b) read-count scaled:    s_rc_i = genome_size / mean_j(C_j), j over the samples in i's pool
 *       ASSUMPTION: every library in a pool started from the same number of cells and
 *       was pooled by equal volume, so differences in fragment yield reflect epitope
 *       abundance. One factor per pool keeps each sample's relative depth:
 *       bw_rc_i = bw_cov_i x C_i / mean_j(C_j). Only meaningful if the assumption holds.
 */

/*
 * Step 9. Fragment BED from the name-sorted BAM (bedtools bamtobed -bedpe), keeping
 * pairs on one chromosome with a fragment < 1000 bp on a canonical chromosome
 * (canonical.chrom.sizes already leaves out chrM when keep_chrM = false).
 * Fragment = leftmost mate start to rightmost mate end.
 * With blacklist_filter = true, fragments overlapping the blacklist are removed.
 * Also writes the raw (unnormalized) fragment-pileup bedGraph that SEACR needs.
 */
process FRAGMENTS {
    tag "${meta.id}${suffix}"
    label 'process_low'
    container params.container__bedtools
    publishDir "${params.outdir}/fragments", mode: 'copy', pattern: '*.fragments.bed.gz'
    publishDir "${params.outdir}/bedgraph/raw", mode: 'copy', pattern: '*.raw.bedGraph'

    input:
    tuple val(meta), path(bam)
    path chrom_sizes
    path blacklist      // assets/NO_FILE when blacklist filtering is off
    val  suffix         // '' for all reads, '.dedup' for the deduplicated control

    output:
    tuple val(meta), path("${meta.id}${suffix}.fragments.bed.gz"), emit: bed
    tuple val(meta), path("${meta.id}${suffix}.raw.bedGraph"),     emit: bedgraph
    tuple val(meta), path("${meta.id}${suffix}.fragments.tsv"),    emit: stats
    path 'versions.yml',                                           emit: versions

    script:
    def out = "${meta.id}${suffix}"
    def filter = blacklist.name != 'NO_FILE'
    """
    export LC_ALL=C
    bedtools bamtobed -bedpe -i ${bam} \\
        | awk -v OFS='\\t' 'NR == FNR { keep[\$1] = 1; next }
                            \$1 == \$4 && (\$1 in keep) {
                                s = \$2 < \$5 ? \$2 : \$5; e = \$3 > \$6 ? \$3 : \$6
                                if (e - s < 1000) print \$1, s, e }' ${chrom_sizes} - \\
        | sort -k1,1 -k2,2n -k3,3n > all.bed

    ${filter ? "bedtools intersect -v -a all.bed -b ${blacklist} > ${out}.fragments.bed" : "mv all.bed ${out}.fragments.bed"}

    bedtools genomecov -bg -i ${out}.fragments.bed -g ${chrom_sizes} > ${out}.raw.bedGraph

    before=\$( [ -f all.bed ] && wc -l < all.bed || wc -l < ${out}.fragments.bed )
    after=\$(wc -l < ${out}.fragments.bed)
    coverage=\$(awk '{ s += \$3 - \$2 } END { printf "%d", s }' ${out}.fragments.bed)
    {
        printf 'fragments_before_blacklist\\t%d\\n' \$before
        printf 'fragments_removed_blacklist\\t%s\\n' "${filter ? '\$(( before - after ))' : 'NA'}"
        printf 'fragments\\t%d\\n' \$after
        printf 'coverage_bp\\t%s\\n' \$coverage
    } > ${out}.fragments.tsv
    rm -f all.bed
    gzip -n ${out}.fragments.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}

// Step 10. scale_factors.tsv: one row per sample, from the fragment stats and genome_size.
process SCALE_FACTORS {
    container params.container__python
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  pools          // list of "<id>=<pool>"
    path stats          // <id>.fragments.tsv of every sample
    path genome_size

    output:
    path 'scale_factors.tsv', emit: tsv

    script:
    """
    scale_factors.py --genome-size ${genome_size} --blacklist-filter ${params.blacklist_filter} \\
        --pools ${pools.join(',')} ${stats} > scale_factors.tsv
    cat scale_factors.tsv
    """
}

// Step 10. Multiply the raw bedGraph by a scale factor and write a bigWig.
// kind = coverage_norm (s_cov) or readcount_scaled (s_rc).
process BIGWIG {
    tag "${meta.id}"
    container params.container__ucsc
    publishDir "${params.outdir}/bigwig/${kind}", mode: 'copy', pattern: '*.bw'

    input:
    tuple val(meta), path(bedgraph), val(factor)
    path chrom_sizes
    val  kind

    output:
    tuple val(meta), path("${meta.id}.${kind}.bw"), emit: bigwig
    path 'versions.yml',                            emit: versions

    script:
    """
    export LC_ALL=C
    awk -v OFS='\\t' -v OFMT='%.7g' -v s=${factor} '{ print \$1, \$2, \$3, \$4 * s }' ${bedgraph} \\
        | sort -k1,1 -k2,2n > scaled.bedGraph
    bedGraphToBigWig scaled.bedGraph ${chrom_sizes} ${meta.id}.${kind}.bw
    rm scaled.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedGraphToBigWig: \$(bedGraphToBigWig 2>&1 | awk '/bedGraphToBigWig v/ { print \$3; exit }')
    END_VERSIONS
    """
}
