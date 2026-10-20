# Workflow steps

This page explains each step of the pipeline in detail. The
[workflow diagram](../README.md#workflow) in the README shows how the steps
connect. Each section below matches one color of the diagram.

- [All steps](#all-steps)
- [What the diagram leaves out](#what-the-diagram-leaves-out)
- [Read preparation](#read-preparation)
- [Alignment and filtering](#alignment-and-filtering)
- [Fragments](#fragments)
- [Coverage and bigWigs](#coverage-and-bigwigs)
- [Feature counts](#feature-counts)
- [Peaks and signal plots](#peaks-and-signal-plots)
- [QC reports](#qc-reports)

## All steps

| Step | Stage | What happens | Tool |
|---|---|---|---|
| 1 | [Read preparation](#read-preparation) | Join the lanes of each sample, then run FastQC | cat, FastQC |
| 2 | [Read preparation](#read-preparation) | Adapter and quality trimming. A pair is dropped if a read becomes shorter than 20 bp | Trim Galore (cutadapt) |
| 3 | [Alignment and filtering](#alignment-and-filtering) | Align to the genome with the Henikoff CUT&Tag settings: `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700` | Bowtie2 |
| 4 | [Alignment and filtering](#alignment-and-filtering) | E. coli carry-over QC: the same reads against E. coli K-12 MG1655 (U00096.3), with `--no-overlap --no-dovetail` added. QC only | Bowtie2 |
| 5 | [Alignment and filtering](#alignment-and-filtering) | Keep properly paired, primary reads with MAPQ >= `min_mapq` (`-f 2 -F 0x904 -q`). Drop a mate whose partner failed the filter. Sort and index | samtools |
| 6 | [Alignment and filtering](#alignment-and-filtering) | Mark duplicates, but keep them (`REMOVE_DUPLICATES=false`) | Picard MarkDuplicates |
| 7 | [QC reports](#qc-reports) | `alignment_stats.tsv`: one row per sample | Python |
| 8 | [QC reports](#qc-reports) | Fragment size counts per sample and one combined plot (0-700 bp) | samtools, matplotlib |
| 9 | [Fragments](#fragments) | Fragment BED: `bedtools bamtobed -bedpe`, same chromosome, fragment < 1000 bp, canonical chromosomes. Optional blacklist filter | bedtools |
| 10 | [Coverage and bigWigs](#coverage-and-bigwigs) | Raw-count bedGraphs (for SEACR), `scale_factors.tsv`, coverage-normalized and (optional) read-count scaled bigWigs | bedtools, bedGraphToBigWig |
| 11 | [Feature counts](#feature-counts) | Count fragments over protein-coding gene bodies, promoters (TSS +/- `promoter_window`) and ENCODE cCREs | featureCounts |
| 12 | [Peaks and signal plots](#peaks-and-signal-plots) | Peaks, stringent mode, with the matched IgG control or the top 1% threshold | SEACR 1.3 |
| 13 | [Peaks and signal plots](#peaks-and-signal-plots) | `peaks_frip.tsv`: peak count and widths, FRiP, control FRiP | bedtools |
| 14 | [Peaks and signal plots](#peaks-and-signal-plots) | Signal at peak summits (per target, with its control) and at protein-coding TSSs (all samples) | deepTools |
| 15 | [Peaks and signal plots](#peaks-and-signal-plots) | Sample correlation: Pearson on log1p values of 5 kb bins, clustered heatmap | deepTools |
| 16 | [QC reports](#qc-reports) | One report of every QC result | MultiQC |

## What the diagram leaves out

The diagram leaves out some inputs, so that it stays readable:

- `QC_TABLES` also reads the logs and stats of steps 1-9.
- `MULTIQC` also reads the trimming, alignment, duplication, featureCounts,
  scale-factor, FRiP and correlation outputs.
- `SOFTWARE_VERSIONS` collects the tool versions of every step into
  `software_versions.yml`.
- `FRAGMENTS` also reads the chromosome sizes and, with
  `blacklist_filter = true`, the blacklist.

## Read preparation

Steps 1 and 2. Several FASTQ pairs per sample (lanes, or a re-sequenced
library) are joined first. FastQC then checks the joined reads. Trim Galore
removes adapters and low-quality bases. A pair is dropped if a read becomes
shorter than 20 bp.

## Alignment and filtering

Steps 3 to 6.

- **Alignment (step 3).** Bowtie2 aligns each pair to the genome with the
  Henikoff CUT&Tag settings:
  `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700`.
- **E. coli carry-over (step 4).** The same reads are aligned to E. coli K-12
  MG1655 (U00096.3), with `--no-overlap --no-dovetail` added. This is QC only.
  The E. coli index is built at each run from `assets/U00096.3.fa.gz`, which
  ships with the pipeline. To skip the build, set `ecoli_index` to a prebuilt
  Bowtie2 index prefix.
- **Filter (step 5).** samtools keeps properly paired, primary reads with
  MAPQ >= `min_mapq` (`-f 2 -F 0x904 -q`). It drops a mate whose partner failed
  the filter. Then it sorts and indexes the BAM.
- **Duplicates (step 6).** Picard MarkDuplicates marks duplicates, but keeps
  them (`REMOVE_DUPLICATES=false`). Duplicates are kept for controls too.

**Optical duplicates** need Illumina read names (with the flow-cell
coordinates). Reads renamed by SRA (for example `SRR12246717.1`) have none, so
the optical duplicate columns are `NA` and the log has a warning. Set
`optical_pixel_distance` to 2500 for patterned flow cells (NovaSeq, HiSeq
X/3000/4000) and to 100 for unpatterned flow cells (HiSeq 2500, NextSeq,
MiSeq).

The genome reference (Bowtie2 index, chromosome sizes) is described in
[Supported genomes](inputs.md#supported-genomes).

## Fragments

Step 9. `bedtools bamtobed -bedpe` turns each read pair into one fragment. A
fragment is kept if both mates are on the same chromosome, it is shorter than
1000 bp and it is on a canonical chromosome. The fragment BED is
`fragments/<sample>.fragments.bed.gz`.

Steps 11-15 all use the fragments from step 9, or the bigWigs made from them.
So the canonical-chromosome restriction and the optional blacklist filter apply
to all of them.

With `dedup_controls = true`, `FRAGMENTS_DEDUP` also makes deduplicated
fragments of each IgG control. Only SEACR uses them. See
[IgG controls](inputs.md#igg-controls).

### Blacklist filtering (optional)

Off by default. With `blacklist_filter = true`:

- Fragments that overlap any ENCODE blacklist v2 region (`bedtools intersect -v`)
  are removed before every later step: bigWigs, SEACR, featureCounts, FRiP,
  deepTools and correlation.
- `genome_size` = canonical genome size - blacklisted bp. A uniform coverage
  still scales to 1.
- `scale_factors.tsv` records the filter state, the fragments removed per
  sample and the `genome_size` used.

With `blacklist_filter = false`, the blacklist file is not read and nothing
about it is reported.

## Coverage and bigWigs

Step 10. This step writes the raw-count bedGraphs (the SEACR input),
`scale_factors.tsv`, the coverage-normalized bigWigs and, optionally, the
read-count scaled bigWigs.

`C_i` is the total per-base coverage of sample i: the sum of its fragment
lengths. `genome_size` is defined in
[Supported genomes](inputs.md#supported-genomes).

**Coverage-normalized bigWigs** (always made, `bigwig/coverage_norm/`):

```
s_cov_i = genome_size / C_i
bw_cov_i(x) = (fragments of sample i covering base x) x s_cov_i
```

Each base holds the fraction of the sample's total coverage at that base, times
`genome_size`. If coverage were uniform, every base would be 1. So the
genome-length-weighted mean of every coverage-normalized bigWig is 1.

**Read-count scaled bigWigs** (only with `readcount_scaled_bw = true`, `bigwig/readcount_scaled/`):

```
s_rc_i  = genome_size / mean_j(C_j)        j = every sample in the same pool
bw_rc_i = bw_cov_i x C_i / mean_j(C_j)
```

The assumption: every library in a pool started from the same number of cells
and was pooled by equal volume. Then a sample with more fragments has more of
the epitope, and this difference should be kept, not normalized away. All
samples of a pool share one scale factor, so each keeps its relative depth. The
average sample of a pool averages about 1 per bp. Controls are included in
their pool.

> **Warning:** the read-count scaled bigWigs are only meaningful if the
> assumption above holds for your libraries. If you are not sure, use the
> coverage-normalized bigWigs.

`scale_factors.tsv` lists, per sample: pool, fragments, `C_i`, `s_cov`,
`s_rc`, `C_i / mean_j(C_j)`, the blacklist filter state, the fragments removed
by the blacklist and the `genome_size` used.

## Feature counts

Step 11. featureCounts counts fragments over protein-coding gene bodies,
promoters (TSS +/- `promoter_window`) and ENCODE cCREs. dm6 has no cCREs, so
cCRE counting is skipped for dm6.

**featureCounts counts fragments, over their full span.** Each fragment from
step 9 becomes one record that covers the whole fragment. A fragment counts for a
feature if it overlaps it by at least 1 bp. A fragment that overlaps several
features counts for each of them (`-O`). This makes the counts use exactly the
same fragments as the bigWigs and FRiP. The spec asked for paired-end
counting of the BAM; counting the fragments instead is a deliberate change.

## Peaks and signal plots

Steps 12 to 15.

- **Peaks (step 12).** SEACR 1.3 calls peaks in stringent mode on targets only.
  A target with a control uses its IgG control; a target without one uses the
  top 1% threshold. The full control rules are in
  [IgG controls](inputs.md#igg-controls).
- **Peak stats (step 13).** `peaks_frip.tsv` has the peak count and widths,
  FRiP and control FRiP of each target.
- **deepTools plots (step 14).** Signal at peak summits (per target, with its
  control) and at protein-coding TSSs (all samples). `plot_window` sets the
  +/- bp around each point.
- **Sample correlation (step 15).** Pearson correlation on log1p values of 5 kb
  bins, with a clustered heatmap.

**Sample correlation uses true 5 kb bin means.** The raw bedGraphs are written
with `bedtools genomecov -bga`, so bases without fragments are explicit 0
intervals in the bigWigs. A deepTools bin value is then the mean over all bases
of the bin. (With `-bg`, deepTools would average over the covered bases only,
which inflates sparse bins several-fold.) Bins that are 0 in every sample are
skipped (`--skipZeros`). The log1p is applied to the bin values before
`plotCorrelation`, because its `--log1p` option changes only scatter plots.

`--skipZeros` keeps only bins with fragments in at least one sample. With
shallow data, most kept bins then have fragments in one sample only, and the
correlations between samples are negative. On the 1M-pair smoke test, every
pair of samples has r < 0. The correlation reflects real co-occupancy only at
higher depth.

## QC reports

Steps 7, 8 and 16.

- **Alignment stats (step 7).** `alignment_stats.tsv` has one row per sample.
  The columns are listed in [Output files](outputs.md#alignment-stats-columns).
- **Fragment sizes (step 8).** Fragment size counts per sample and one combined
  plot (0-700 bp).
- **MultiQC (step 16).** One report of every QC result.
