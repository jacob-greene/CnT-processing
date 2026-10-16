/*
 * Steps 3-6: Bowtie2 alignment (target genome and E. coli), read-pair filtering,
 * fragment-length counts and duplicate marking.
 */

// Henikoff lab CUT&Tag settings (Zheng et al., CUT&Tag data processing tutorial).
def bt2CutTag() { '--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700' }

// Step 3. Align to the target genome. Output: unsorted BAM + the Bowtie2 log.
process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_medium'
    container params.container__bowtie2
    publishDir "${params.outdir}/alignment/bowtie2_logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(r1), path(r2)
    tuple path(index_dir), val(index_name)

    output:
    tuple val(meta), path("${meta.id}.bam"),         emit: bam
    tuple val(meta), path("${meta.id}.bowtie2.log"), emit: log
    path 'versions.yml',                             emit: versions

    script:
    """
    bowtie2 ${bt2CutTag()} -p ${task.cpus} -x ${index_dir}/${index_name} \\
            -1 ${r1} -2 ${r2} 2> ${meta.id}.bowtie2.log \\
        | samtools view -b -o ${meta.id}.bam -
    cat ${meta.id}.bowtie2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | awk 'NR == 1 { print \$NF }')
        samtools: \$(samtools --version | head -1 | awk '{ print \$2 }')
    END_VERSIONS
    """
}

// Step 4. E. coli carry-over QC: same reads, same settings plus --no-overlap --no-dovetail.
// Only the log is kept: the concordant pair count is the carry-over. No spike-in normalization.
process BOWTIE2_ECOLI {
    tag "${meta.id}"
    label 'process_medium'
    container params.container__bowtie2
    publishDir "${params.outdir}/alignment/bowtie2_logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(r1), path(r2)
    tuple path(index_dir), val(index_name)

    output:
    tuple val(meta), path("${meta.id}.ecoli.bowtie2.log"), emit: log

    script:
    """
    bowtie2 ${bt2CutTag()} --no-overlap --no-dovetail -p ${task.cpus} -x ${index_dir}/${index_name} \\
            -1 ${r1} -2 ${r2} -S /dev/null 2> ${meta.id}.ecoli.bowtie2.log
    cat ${meta.id}.ecoli.bowtie2.log
    """
}

/*
 * Step 5. Keep properly paired, primary alignments with MAPQ >= min_mapq
 * (samtools view -f 2 -F 0x904 -q min_mapq), then sort and index.
 *
 * Bowtie2 gives each mate its own MAPQ, so the MAPQ filter can keep one mate and
 * drop the other. Orphan mates are removed after the filter, so every kept read
 * has its mate and "filtered pairs" = kept reads / 2.
 *
 * Also written here, from the filtered BAM (before any chrM exclusion or blacklist
 * filtering): samtools flagstat and idxstats, the pairs on mito_contig, the
 * fragment length counts (|TLEN| of read 1), and whether read names hold
 * Illumina flow-cell coordinates (needed for optical duplicate detection).
 */
process FILTER_BAM {
    tag "${meta.id}"
    label 'process_low'
    container params.container__samtools
    publishDir "${params.outdir}/alignment/stats", mode: 'copy', pattern: '*.{flagstat,idxstats}'
    publishDir "${params.outdir}/fragment_sizes", mode: 'copy', pattern: '*.fragment_sizes.tsv',
        saveAs: { "${meta.id}.tsv" }

    input:
    tuple val(meta), path(bam)
    val mito_contig

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.filter.tsv"),             emit: stats
    tuple val(meta), path("${meta.id}.fragment_sizes.tsv"),     emit: fragment_sizes
    tuple val(meta), path("${meta.id}.filtered.flagstat"),      emit: flagstat
    tuple val(meta), path("${meta.id}.filtered.idxstats"),      emit: idxstats
    path 'versions.yml',                                        emit: versions

    script:
    """
    samtools view -u -f 2 -F 0x904 -q ${params.min_mapq} ${bam} \\
        | samtools sort -n -@ ${task.cpus} -T tmp_n -O sam - \\
        | awk 'BEGIN { FS = "\\t" }
               /^@/ { print; next }
               # Name-sorted: a kept pair is two consecutive records with the same name.
               \$1 == prev { print held; print; prev = ""; next }
               { held = \$0; prev = \$1 }' \\
        | samtools sort -@ ${task.cpus} -T tmp_c -o ${meta.id}.filtered.bam -
    samtools index ${meta.id}.filtered.bam
    samtools flagstat ${meta.id}.filtered.bam > ${meta.id}.filtered.flagstat
    samtools idxstats ${meta.id}.filtered.bam > ${meta.id}.filtered.idxstats

    # Fragment length = |TLEN| of read 1 (flag 64) of each kept pair.
    printf 'length\\tcount\\n' > ${meta.id}.fragment_sizes.tsv
    samtools view -f 64 ${meta.id}.filtered.bam \\
        | awk '{ t = \$9 < 0 ? -\$9 : \$9; n[t]++ } END { for (t in n) print t "\\t" n[t] }' \\
        | sort -k1,1n >> ${meta.id}.fragment_sizes.tsv

    reads=\$(samtools view -c ${meta.id}.filtered.bam)
    mito=\$(awk -v m=${mito_contig} '\$1 == m { print \$3 }' ${meta.id}.filtered.idxstats)
    # Picard reads tile:x:y from the last three ':' fields of a 5- or 7-field Illumina name.
    # SRA-renamed reads (e.g. "SRR12246717.1") have none, so optical duplicates cannot be found.
    first=\$(samtools view ${meta.id}.filtered.bam | head -n 1 | cut -f 1 || true)
    illumina=\$(echo "\$first" | awk -F ':' '(NF == 5 || NF == 7) && \$(NF-2) ~ /^[0-9]+\$/ && \$(NF-1) ~ /^[0-9]+\$/ && \$NF ~ /^[0-9]+\$/ { print "true"; exit } { print "false" }')
    {
        printf 'filtered_pairs\\t%s\\n' \$(( reads / 2 ))
        printf 'mito_pairs\\t%s\\n' \$(( \${mito:-0} / 2 ))
        printf 'first_read_name\\t%s\\n' "\$first"
        printf 'illumina_read_names\\t%s\\n' "\$illumina"
    } > ${meta.id}.filter.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{ print \$2 }')
    END_VERSIONS
    """
}

/*
 * Step 6. Picard MarkDuplicates, duplicates KEPT (REMOVE_DUPLICATES=false).
 * Downstream steps use all reads. If read names have no Illumina coordinates,
 * optical duplicate detection is switched off (READ_NAME_REGEX=null) and the
 * optical columns are reported as NA, not 0.
 */
process MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_low'
    container params.container__picard
    publishDir "${params.outdir}/duplicates", mode: 'copy', pattern: '*.markdup_metrics.txt'

    input:
    tuple val(meta), path(bam), path(bai), val(illumina_names)

    output:
    tuple val(meta), path("${meta.id}.markdup.bam"), emit: bam
    tuple val(meta), path("${meta.id}.markdup_metrics.txt"), val(illumina_names), emit: metrics
    path 'versions.yml', emit: versions

    script:
    def xmx = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    def regex = illumina_names ? '' : 'READ_NAME_REGEX=null'
    """
    # The biocontainer has a `picard` wrapper; the Lmod module only ships picard.jar.
    if command -v picard > /dev/null; then
        run_picard() { picard -Xmx${xmx}g "\$@"; }
    else
        run_picard() { java -Xmx${xmx}g -jar "\$EBROOTPICARD/picard.jar" "\$@"; }
    fi
    run_picard MarkDuplicates \\
        I=${bam} O=${meta.id}.markdup.bam M=${meta.id}.markdup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=${params.optical_pixel_distance} \\
        ${regex} VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(run_picard MarkDuplicates --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}

// Name-sort the marked BAM for bedtools bamtobed -bedpe. With dedup_controls = true,
// controls also get a name-sorted BAM without duplicates (flag 1024), used only for SEACR.
process NAMESORT {
    tag "${meta.id}"
    label 'process_low'
    container params.container__samtools

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.namesorted.bam"),       emit: bam
    tuple val(meta), path("${meta.id}.dedup.namesorted.bam"), emit: dedup, optional: true

    script:
    def dedup = params.dedup_controls && meta.is_control
    """
    samtools sort -n -@ ${task.cpus} -T tmp -o ${meta.id}.namesorted.bam ${bam}
    ${dedup ? "samtools view -b -F 1024 -o ${meta.id}.dedup.namesorted.bam ${meta.id}.namesorted.bam" : ''}
    """
}
