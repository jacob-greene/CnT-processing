/*
 * Steps 1-2: concatenate lanes, FastQC, adapter and quality trimming.
 */

// Step 1a. Join the lanes of one sample into one R1 and one R2 file, and count the read pairs.
process CAT_FASTQ {
    tag "${meta.id}"
    container params.container__samtools

    input:
    tuple val(meta), path(r1, stageAs: 'in_R1/*'), path(r2, stageAs: 'in_R2/*')

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.read_pairs.tsv"),                              emit: stats

    script:
    // Gzip files can be joined with cat. List order = sample sheet order, the same for R1 and R2.
    def files1 = (r1 instanceof List ? r1 : [r1]).join(' ')
    def files2 = (r2 instanceof List ? r2 : [r2]).join(' ')
    """
    cat ${files1} > ${meta.id}_R1.fastq.gz
    cat ${files2} > ${meta.id}_R2.fastq.gz

    n1=\$(( \$(zcat ${meta.id}_R1.fastq.gz | wc -l) / 4 ))
    n2=\$(( \$(zcat ${meta.id}_R2.fastq.gz | wc -l) / 4 ))
    if [ "\$n1" -ne "\$n2" ]; then
        echo "ERROR: ${meta.id}: R1 has \$n1 reads but R2 has \$n2 reads" >&2
        exit 1
    fi
    printf 'total_pairs\\t%s\\n' "\$n1" > ${meta.id}.read_pairs.tsv
    """
}

// Step 1b. FastQC on the joined, untrimmed reads.
process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    container params.container__fastqc
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: '*_fastqc.*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path '*_fastqc.{zip,html}', emit: reports
    path 'versions.yml',        emit: versions

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}

// Step 2. Trim Galore (cutadapt): adapter and quality trimming in paired mode.
// Pairs where either read becomes shorter than 20 bp are dropped (Trim Galore default).
// For the 25 x 25 bp test reads this removes almost nothing.
process TRIM_GALORE {
    tag "${meta.id}"
    label 'process_low'
    container params.container__trimgalore
    publishDir "${params.outdir}/trimming", mode: 'copy', pattern: '*_trimming_report.txt'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_val_1.fq.gz"), path("${meta.id}_val_2.fq.gz"), emit: reads
    path '*_trimming_report.txt', emit: reports
    path 'versions.yml',          emit: versions

    script:
    """
    trim_galore --paired --gzip --basename ${meta.id} ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(trim_galore --version | awk '/version/ { print \$NF; exit }')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
