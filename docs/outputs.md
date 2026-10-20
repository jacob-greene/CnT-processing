# Output files

- [Output folder](#output-folder)
- [Alignment stats columns](#alignment-stats-columns)

Back to the [README](../README.md). How each file is made:
[Workflow steps](steps.md).

## Output folder

```
alignment_stats.tsv            One row per sample: read pairs, alignment rate, filtered pairs,
                               % chrM, duplication %, optical duplicate %, E. coli %, final fragments
duplication.tsv                Picard results: pairs, PERCENT_DUPLICATION, optical duplicates,
                               ESTIMATED_LIBRARY_SIZE
scale_factors.tsv              C_i, s_cov, s_rc and the genome_size used, per sample
peaks_frip.tsv                 One row per target: peaks, widths, FRiP, control_frip, top0.01 set
bigwig/coverage_norm/          <sample>.coverage_norm.bw
bigwig/readcount_scaled/       <sample>.readcount_scaled.bw (readcount_scaled_bw = true)
bedgraph/raw/                  <sample>.raw.bedGraph: unnormalized fragment pileup (SEACR input)
fragments/                     <sample>.fragments.bed.gz: the fragments every later step uses
fragment_sizes/                <sample>.tsv (length, count) and fragment_sizes.pdf / .png
featurecounts/                 genes / promoters / ccre .counts.tsv (features x samples) and .summary
featurecounts/annotation/      The SAF files that were counted, and tss_genes.bed
seacr/                         <target>.stringent.bed: primary peaks
seacr/top0.01/                 <target>.stringent.bed: numeric-mode peaks (targets with a control)
deeptools/peaks/               <target>.peaks_heatmap.pdf, .peaks_profile.pdf (with the control)
deeptools/tss/                 tss_heatmap_coverage_norm.pdf, tss_profile_coverage_norm.pdf,
                               tss_profile_readcount_scaled.pdf (same y-axis range)
correlation/                   correlation_heatmap.pdf / .png, correlation_matrix.tsv,
                               correlation_bins.tsv.gz (5 kb bin values)
multiqc/multiqc_report.html    The MultiQC report
alignment/                     Bowtie2 logs (target and E. coli), samtools flagstat and idxstats
duplicates/                    Picard MarkDuplicates metrics
fastqc/, trimming/             FastQC reports and Trim Galore reports
genome/                        canonical.chrom.sizes, genome_size.tsv, genes.gtf.gz
pipeline_info/                 software_versions.yml
```

## Alignment stats columns

`alignment_stats.tsv` columns:

| Column | Meaning |
|---|---|
| `total_pairs` | Read pairs in the joined FASTQs |
| `pairs_after_trimming` | Read pairs after Trim Galore |
| `aligned_pairs` | Pairs that Bowtie2 aligned concordantly |
| `alignment_rate_pct` | `aligned_pairs / pairs_after_trimming` x 100 (the pairs Bowtie2 was given) |
| `filtered_pairs` | Pairs that pass the MAPQ, proper-pair and primary filter |
| `pct_chrM` | Filtered pairs on `chrM` / all filtered pairs x 100. Computed before any chrM or blacklist exclusion |
| `duplication_pct` | Picard `PERCENT_DUPLICATION` x 100 (duplicates are kept) |
| `optical_duplicate_pct` | Optical duplicate pairs / filtered pairs x 100, or `NA` |
| `ecoli_pairs`, `ecoli_pct` | Pairs that align concordantly to E. coli, and their % of the pairs after trimming |
| `final_fragments` | Fragments used by every later step |
