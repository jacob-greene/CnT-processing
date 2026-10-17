/*
 * Steps 6-8 summary tables and plots: alignment_stats.tsv, duplication.tsv,
 * the combined fragment size plot, and software_versions.yml.
 */

// Steps 6-7. One row per sample. qc_tables.py also checks that the pair counts reconcile.
process QC_TABLES {
    container params.container__python
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  ids            // sample ids, in sample sheet order
    path files          // every per-sample stats file and log, named <id>.<kind>

    output:
    path 'alignment_stats.tsv', emit: alignment_stats
    path 'duplication.tsv',     emit: duplication

    script:
    """
    qc_tables.py ${ids.join(' ')}
    column -t alignment_stats.tsv 2> /dev/null || cat alignment_stats.tsv
    """
}

// Step 8. All samples' fragment size distributions on one plot (0-700 bp).
process FRAGMENT_SIZE_PLOT {
    container params.container__python
    publishDir "${params.outdir}/fragment_sizes", mode: 'copy'

    input:
    path tables

    output:
    path 'fragment_sizes.{pdf,png}', emit: plot

    script:
    """
    plot_fragment_sizes.py ${tables}
    """
}

// Merge the versions.yml of every process into software_versions.yml.
process SOFTWARE_VERSIONS {
    container params.container__python
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path versions, stageAs: 'versions/v*.yml'

    output:
    path 'software_versions.yml'

    script:
    """
    {
        echo 'Workflow:'
        echo '    CnT-processing: ${workflow.manifest.version}'
        echo '    Nextflow: ${nextflow.version}'
        # One block per process; identical blocks (one per sample) are kept once.
        cat versions/*.yml | awk '/^"/ { key = \$0; first = !(key in blk); if (first) blk[key] = ""; next }
                                  first { blk[key] = blk[key] \$0 "\\n" }
                                  END   { for (k in blk) printf "%s\\n%s", k, blk[k] }'
    } > software_versions.yml
    cat software_versions.yml
    """
}

/*
 * Step 16. MultiQC over FastQC, Trim Galore, Bowtie2 (target and E. coli), samtools
 * flagstat, Picard duplicates and featureCounts, plus the pipeline's own tables
 * (alignment_stats, peaks_frip, duplication, scale_factors), the fragment sizes and
 * the correlation heatmap as custom content. Module settings: assets/multiqc_config.yml.
 */
process MULTIQC {
    label 'process_low'
    container params.container__multiqc
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path config
    path 'inputs/*'
    path 'tables/*'
    path 'fragment_sizes/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data',        emit: data
    path 'versions.yml',        emit: versions

    script:
    """
    # The pipeline's own tables as MultiQC custom-content tables.
    for t in alignment_stats duplication scale_factors peaks_frip; do
        [ -f tables/\$t.tsv ] || continue
        title=\$(echo \$t | tr '_' ' ')
        { printf "# id: '%s'\\n# section_name: '%s'\\n# plot_type: 'table'\\n" "\$t" "\$title"
          cat tables/\$t.tsv; } > inputs/\${t}_mqc.tsv
    done
    [ -f tables/correlation_heatmap.png ] && cp tables/correlation_heatmap.png inputs/correlation_heatmap_mqc.png
    fragment_sizes_mqc.py fragment_sizes/*.tsv > inputs/fragment_sizes_mqc.json

    multiqc --force --config ${config} --title 'CnT-processing' --filename multiqc_report.html inputs
    # --filename also renames the data folder (multiqc_report_data); keep the usual name.
    [ -d multiqc_report_data ] && mv multiqc_report_data multiqc_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.* version //')
    END_VERSIONS
    """
}
