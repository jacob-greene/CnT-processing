/*
 * Steps 11-15: featureCounts, SEACR peaks, peak stats + FRiP, deepTools plots
 * and sample correlation.
 *
 * Every step reads the fragment files from step 9 (or the bigWigs made from
 * them), so the canonical-chromosome restriction and the optional blacklist
 * filter apply to all of them.
 */

// Step 11 input. Protein-coding gene bodies, promoters (TSS +/- promoter_window)
// and cCREs as SAF, plus the gene BED6 used for the TSS plots.
process ANNOTATIONS {
    tag "${params.genome}"
    container params.container__python
    publishDir "${params.outdir}/featurecounts/annotation", mode: 'copy'

    input:
    path gtf
    path chrom_sizes
    path ccre           // assets/NO_FILE when the genome has no cCREs (dm6)

    output:
    path '*.saf',                  emit: saf
    path 'tss_genes.bed',          emit: tss_bed
    path 'annotation_summary.tsv', emit: summary

    script:
    def ccre_arg = ccre.name != 'NO_FILE' ? "--ccre ${ccre}" : ''
    """
    make_annotations.py --gtf ${gtf} --chrom-sizes ${chrom_sizes} \\
        --promoter-window ${params.promoter_window} ${ccre_arg}
    cat annotation_summary.tsv
    """
}

/*
 * Step 11. featureCounts: one run per feature set over all samples, giving one
 * matrix (features x samples, controls included) and one summary.
 *
 * Each fragment from the fragment BED is counted once, over its full span: the
 * BED is turned into a SAM file with one record per fragment (CIGAR = fragment
 * length). This counts the same fragments as the bigWigs and FRiP (canonical
 * chromosomes, < 1000 bp, blacklist-filtered when blacklist_filter = true).
 * -O: a fragment that overlaps several features is counted for each of them.
 */
process FEATURECOUNTS {
    tag "${set_name}"
    label 'process_low'
    container params.container__subread
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    tuple val(set_name), path(saf)
    val  ids            // sample ids, in sample sheet order (= matrix column order)
    path beds           // <id>.fragments.bed.gz of every sample
    path chrom_sizes

    output:
    path "${set_name}.counts.tsv",         emit: counts
    path "${set_name}.counts.tsv.summary", emit: summary
    path 'versions.yml',                   emit: versions

    script:
    """
    for id in ${ids.join(' ')}; do
        zcat \$id.fragments.bed.gz \\
            | awk -v OFS='\\t' 'NR == FNR { print "@SQ", "SN:" \$1, "LN:" \$2; next }
                                { print "f" FNR, 0, \$1, \$2 + 1, 255, (\$3 - \$2) "M", "*", 0, 0, "*", "*" }' \\
                  ${chrom_sizes} - > \$id.sam
    done
    featureCounts -F SAF -a ${saf} -O -T ${task.cpus} -o raw.tsv \$(for id in ${ids.join(' ')}; do echo \$id.sam; done)

    # Column names: "<id>.sam" -> "<id>"; drop the command-line comment line.
    grep -v '^#' raw.tsv | sed '1 s/\\.sam//g' > ${set_name}.counts.tsv
    sed '1 s/\\.sam//g' raw.tsv.summary > ${set_name}.counts.tsv.summary
    rm -f *.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+' | sed 's/^v//')
    END_VERSIONS
    """
}

/*
 * Step 12. SEACR 1.3 on the raw fragment-pileup bedGraphs (stringent mode).
 *   target with a control:    SEACR <target.bg> <control.bg> norm stringent  -> primary peaks
 *                             SEACR <target.bg> 0.01 non stringent           -> seacr/top0.01/
 *   target without a control: SEACR <target.bg> 0.01 non stringent           -> primary peaks
 * With dedup_controls = true, the control bedGraph is the deduplicated one.
 * Controls are not peak-called.
 */
process SEACR {
    tag "${meta.id}"
    container params.container__seacr
    publishDir "${params.outdir}/seacr", mode: 'copy', pattern: "${meta.id}.stringent.bed"
    publishDir "${params.outdir}/seacr/top0.01", mode: 'copy', pattern: "top0.01/*.bed",
        saveAs: { new File(it).name }

    input:
    tuple val(meta), path(bedgraph), path(control_bedgraph)    // control: assets/NO_FILE if none

    output:
    tuple val(meta), path("${meta.id}.stringent.bed"), path("top0.01/*"), val(mode), emit: peaks
    path 'versions.yml', emit: versions

    script:
    def has_control = control_bedgraph.name != 'NO_FILE'
    mode = has_control ? 'control' : 'top0.01'
    """
    # SEACR is not an Lmod module: outside the container, use the pinned v1.3 scripts from GitHub.
    if command -v SEACR_1.3.sh > /dev/null; then
        seacr() { SEACR_1.3.sh "\$@"; }
    else
        url=https://raw.githubusercontent.com/FredHutch/SEACR/5a0efe59f06fb17cf9d34d415bb0c1a1f7a77a3c
        mkdir -p seacr_bin
        curl -sfL --retry 3 -o seacr_bin/SEACR_1.3.sh \$url/SEACR_1.3.sh
        curl -sfL --retry 3 -o seacr_bin/SEACR_1.3.R  \$url/SEACR_1.3.R
        seacr() { bash seacr_bin/SEACR_1.3.sh "\$@"; }
    fi

    mkdir -p top0.01
    if [ "${has_control}" = true ]; then
        seacr ${bedgraph} ${control_bedgraph} norm stringent ${meta.id}
        seacr ${bedgraph} 0.01 non stringent top0.01/${meta.id}
        touch top0.01/${meta.id}.stringent.bed      # SEACR writes no file when it finds no peak
    else
        seacr ${bedgraph} 0.01 non stringent ${meta.id}
        touch top0.01/NONE                          # the primary set is already the 0.01 set
    fi
    touch ${meta.id}.stringent.bed
    rm -rf seacr_bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seacr: 1.3
    END_VERSIONS
    """
}

/*
 * Step 13. One peaks_frip.tsv row per target, from the primary peaks and (when a
 * control exists) the top0.01 peaks:
 *   FRiP         = target fragments overlapping >= 1 peak / all target fragments
 *   control_frip = control fragments overlapping >= 1 of the TARGET's peaks /
 *                  all control fragments
 */
process PEAK_STATS {
    tag "${meta.id}"
    container params.container__bedtools

    input:
    tuple val(meta), path(peaks), path(top_peaks, stageAs: 'top/*'), val(mode), path(fragments), path(control_fragments)

    output:
    path "${meta.id}.peaks_frip.tsv", emit: tsv

    script:
    def ctrl = control_fragments.name != 'NO_FILE' ? control_fragments : ''
    def top = top_peaks.name != 'NONE' ? top_peaks : ''
    """
    peak_stats.sh ${meta.id} ${meta.control ?: 'NA'} ${mode} ${fragments} ${peaks} "${ctrl}" "${top}" > ${meta.id}.peaks_frip.tsv
    cat ${meta.id}.peaks_frip.tsv
    """
}

/*
 * Step 14a. Signal at peaks: computeMatrix reference-point --referencePoint center
 * on the SEACR summit regions (column 6), +/- plot_window, on the coverage-normalized
 * bigWigs of the target and its matched control. A target with no peaks gets no plot.
 */
process DEEPTOOLS_PEAKS {
    tag "${meta.id}"
    label 'process_low'
    container params.container__python
    publishDir "${params.outdir}/deeptools/peaks", mode: 'copy'

    input:
    tuple val(meta), path(peaks), path(bigwig), path(control_bigwig)

    output:
    path "${meta.id}.peaks_*", optional: true, emit: plots
    path 'versions.yml',                       emit: versions

    script:
    def ctrl = control_bigwig.name != 'NO_FILE'
    def labels = ctrl ? "${meta.id} ${meta.control}" : "${meta.id}"
    """
    # Column 6 of a SEACR peak = the region with the maximum signal ("chr:start-end").
    awk -v OFS='\\t' '{ split(\$6, a, /[:-]/); print a[1], a[2], a[3], "peak" NR }' ${peaks} > summits.bed
    if [ -s summits.bed ]; then
        computeMatrix reference-point --referencePoint center -R summits.bed \\
            -S ${bigwig} ${ctrl ? control_bigwig : ''} --samplesLabel ${labels} \\
            -a ${params.plot_window} -b ${params.plot_window} --missingDataAsZero \\
            -p ${task.cpus} -o ${meta.id}.peaks_matrix.gz
        plotHeatmap -m ${meta.id}.peaks_matrix.gz -o ${meta.id}.peaks_heatmap.pdf \\
            --refPointLabel summit --regionsLabel "${meta.id} peaks"
        plotProfile -m ${meta.id}.peaks_matrix.gz -o ${meta.id}.peaks_profile.pdf --perGroup \\
            --refPointLabel summit --regionsLabel "${meta.id} peaks"
        plotHeatmap -m ${meta.id}.peaks_matrix.gz -o ${meta.id}.peaks_heatmap.png \\
            --refPointLabel summit --regionsLabel "${meta.id} peaks"
    else
        echo "${meta.id}: no SEACR peaks, so no peak plots." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed 's/computeMatrix //')
    END_VERSIONS
    """
}

/*
 * Step 14b. Signal at TSSs: computeMatrix reference-point --referencePoint TSS on the
 * protein-coding genes (strand-aware), +/- plot_window, all samples in one heatmap and
 * one profile. With readcount_scaled_bw = true, the profile is repeated on the
 * read-count scaled bigWigs, and both profiles share the same y-axis limits.
 */
process DEEPTOOLS_TSS {
    label 'process_low'
    container params.container__python
    publishDir "${params.outdir}/deeptools/tss", mode: 'copy'

    input:
    val  ids            // sample ids, in sample sheet order
    path cov_bigwigs, stageAs: 'coverage_norm/*'
    path rc_bigwigs,  stageAs: 'readcount_scaled/*'     // empty list when readcount_scaled_bw = false
    path tss_bed

    output:
    path 'tss_*', emit: plots
    path 'versions.yml', emit: versions

    script:
    def kinds = rc_bigwigs ? ['coverage_norm', 'readcount_scaled'] : ['coverage_norm']
    """
    for kind in ${kinds.join(' ')}; do
        computeMatrix reference-point --referencePoint TSS -R ${tss_bed} \\
            -S \$(for id in ${ids.join(' ')}; do echo \$kind/\$id.\$kind.bw; done) \\
            --samplesLabel ${ids.join(' ')} -a ${params.plot_window} -b ${params.plot_window} \\
            --missingDataAsZero -p ${task.cpus} -o tss_matrix_\$kind.gz
        plotProfile -m tss_matrix_\$kind.gz -o profile.png --outFileNameData tss_profile_\$kind.tsv
    done
    # One y-axis range for every profile (min and max of every mean signal value).
    range=\$(cat tss_profile_*.tsv | awk -F '\\t' '\$1 != "bin labels" && \$1 != "bins" {
                for (i = 3; i <= NF; i++) if (\$i != "") { v = \$i + 0; if (!n++ || v < lo) lo = v; if (n == 1 || v > hi) hi = v } }
            END { pad = (hi - lo) * 0.05; printf "%g %g", lo - pad, hi + pad }')
    set -- \$range
    for kind in ${kinds.join(' ')}; do
        for ext in pdf png; do
            plotProfile -m tss_matrix_\$kind.gz -o tss_profile_\$kind.\$ext --perGroup \\
                --yMin \$1 --yMax \$2 --plotTitle "TSS profile (\$kind)" --regionsLabel "protein-coding genes"
        done
    done
    plotHeatmap -m tss_matrix_coverage_norm.gz -o tss_heatmap_coverage_norm.pdf --regionsLabel "protein-coding genes"
    plotHeatmap -m tss_matrix_coverage_norm.gz -o tss_heatmap_coverage_norm.png --regionsLabel "protein-coding genes"
    rm -f profile.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed 's/computeMatrix //')
    END_VERSIONS
    """
}

/*
 * Step 15. Sample correlation: multiBigwigSummary bins (5 kb) on the coverage-normalized
 * bigWigs, then plotCorrelation: Pearson on log1p values, hierarchically clustered heatmap.
 * Controls included.
 *
 * A bin without fragments has no bigWig value (NaN). plotCorrelation would drop every bin
 * that is NaN in ANY sample, which keeps only bins covered in every sample (with a sparse
 * IgG, mostly chrM). So NaN is set to 0 (no fragments = zero coverage), and bins that are
 * 0 in every sample are skipped (--skipZeros).
 *
 * plotCorrelation --log1p changes only the scatter plot axes, not the correlation. So the
 * log1p is applied to the bin matrix here, before plotCorrelation reads it.
 */
process CORRELATION {
    label 'process_low'
    container params.container__python
    publishDir "${params.outdir}/correlation", mode: 'copy'

    input:
    val  ids
    path bigwigs

    output:
    path 'correlation_matrix.tsv',      emit: matrix
    path 'correlation_heatmap.{pdf,png}', emit: plots
    path 'correlation_bins.tsv.gz',     emit: bins
    path 'versions.yml',                emit: versions

    script:
    """
    multiBigwigSummary bins -b \$(for id in ${ids.join(' ')}; do echo \$id.coverage_norm.bw; done) \\
        --labels ${ids.join(' ')} --binSize 5000 -p ${task.cpus} \\
        -o bins_nan.npz --outRawCounts correlation_bins.tsv
    python3 -c "import numpy as np; d = dict(np.load('bins_nan.npz', allow_pickle=True)); d['matrix'] = np.log1p(np.nan_to_num(d['matrix'])); np.savez_compressed('bins.npz', **d)"
    for ext in pdf png; do
        plotCorrelation -in bins.npz --corMethod pearson --skipZeros --whatToPlot heatmap \\
            --plotNumbers --colorMap RdYlBu_r --plotTitle "Pearson correlation, log1p, 5 kb bins" \\
            -o correlation_heatmap.\$ext --outFileCorMatrix correlation_matrix.tsv
    done
    gzip -n correlation_bins.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed 's/multiBigwigSummary //')
    END_VERSIONS
    """
}
