# CnT-processing

A small, readable Nextflow (DSL2) pipeline for paired-end CUT&Tag data.
It takes FASTQ files and produces:

1. Coverage-normalized bigWigs
2. Optional read-count scaled bigWigs, for comparison with (1)
3. featureCounts matrices over protein-coding genes, promoters and ENCODE cCREs
4. SEACR peaks (IgG-controlled when a control is given), FRiP, peak stats and deepTools plots at peaks and TSSs
5. QC tables: alignment stats (with % chrM), fragment sizes, duplication (with optical duplicates), E. coli carry-over, sample correlation
6. A MultiQC report

The tool choices and settings follow the Henikoff lab
[CUT&Tag tutorial](https://yezhengstat.github.io/CUTTag_tutorial/) and
[nf-core/cutandrun 1.1](https://nf-co.re/cutandrun/1.1/). This pipeline is not
a fork of nf-core/cutandrun. It does no spike-in normalization.

- [Workflow steps](#workflow-steps)
- [Input data](#input-data)
- [Sample sheet](#sample-sheet)
- [IgG controls](#igg-controls)
- [Supported genomes](#supported-genomes)
- [Normalization](#normalization)
- [Blacklist filtering (optional)](#blacklist-filtering-optional)
- [Parameters](#parameters)
- [Launching with `nextflow run`](#launching-with-nextflow-run)
- [Quickstart with the BASH Workbench](#quickstart-with-the-bash-workbench)
- [Resuming and cleaning up](#resuming-and-cleaning-up)
- [Output files](#output-files)
- [Testing](#testing)
- [Authors](#authors)

## Workflow steps

| Step | What happens | Tool |
|---|---|---|
| 1 | Join the lanes of each sample, then run FastQC | cat, FastQC |
| 2 | Adapter and quality trimming. A pair is dropped if a read becomes shorter than 20 bp | Trim Galore (cutadapt) |
| 3 | Align to the genome with the Henikoff CUT&Tag settings: `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700` | Bowtie2 |
| 4 | E. coli carry-over QC: the same reads against E. coli K-12 MG1655 (U00096.3), with `--no-overlap --no-dovetail` added. QC only | Bowtie2 |
| 5 | Keep properly paired, primary reads with MAPQ >= `min_mapq` (`-f 2 -F 0x904 -q`). Drop a mate whose partner failed the filter. Sort and index | samtools |
| 6 | Mark duplicates, but keep them (`REMOVE_DUPLICATES=false`) | Picard MarkDuplicates |
| 7 | `alignment_stats.tsv`: one row per sample | Python |
| 8 | Fragment size counts per sample and one combined plot (0-700 bp) | samtools, matplotlib |
| 9 | Fragment BED: `bedtools bamtobed -bedpe`, same chromosome, fragment < 1000 bp, canonical chromosomes. Optional blacklist filter | bedtools |
| 10 | Raw-count bedGraphs (for SEACR), `scale_factors.tsv`, coverage-normalized and (optional) read-count scaled bigWigs | bedtools, bedGraphToBigWig |
| 11 | Count fragments over protein-coding gene bodies, promoters (TSS +/- `promoter_window`) and ENCODE cCREs | featureCounts |
| 12 | Peaks, stringent mode, with the matched IgG control or the top 1% threshold | SEACR 1.3 |
| 13 | `peaks_frip.tsv`: peak count and widths, FRiP, control FRiP | bedtools |
| 14 | Signal at peak summits (per target, with its control) and at protein-coding TSSs (all samples) | deepTools |
| 15 | Sample correlation: Pearson on log1p values of 5 kb bins, clustered heatmap | deepTools |
| 16 | One report of every QC result | MultiQC |

Steps 11-15 all use the fragments from step 9, or the bigWigs made from them.
So the canonical-chromosome restriction and the optional blacklist filter apply
to all of them.

**featureCounts counts fragments, over their full span.** Each fragment from
step 9 becomes one record that covers the whole fragment. A fragment counts for a
feature if it overlaps it by at least 1 bp. A fragment that overlaps several
features counts for each of them (`-O`). This makes the counts use exactly the
same fragments as the bigWigs and FRiP. The spec asked for paired-end
counting of the BAM; counting the fragments instead is a deliberate change.

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

## Input data

Paired-end CUT&Tag FASTQ files (gzipped), any read length. Several FASTQ pairs
per sample (lanes, or a re-sequenced library) are joined before any other step.

## Sample sheet

A CSV file with one row per FASTQ pair:

```
sample,fastq_1,fastq_2,control,is_control,pool
H3K27me3_rep1,/path/to/H3K27me3_rep1_R1.fastq.gz,/path/to/H3K27me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1
H3K4me3_rep1,/path/to/H3K4me3_rep1_R1.fastq.gz,/path/to/H3K4me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1
IgG_rep1,/path/to/run1/IgG_rep1_R1.fastq.gz,/path/to/run1/IgG_rep1_R2.fastq.gz,,true,pool1
IgG_rep1,/path/to/run2/IgG_rep1_R1.fastq.gz,/path/to/run2/IgG_rep1_R2.fastq.gz,,true,pool1
```

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample name: letters, digits, `.`, `_` and `-`. Rows with the same name are lanes and are joined, in sheet order. |
| `fastq_1`, `fastq_2` | yes | R1 and R2 FASTQ. A relative path is relative to the folder of the sample sheet. |
| `control` | no | Name of an IgG (negative-control) sample in the same sheet. |
| `is_control` | no | `true` for IgG / negative controls. Default `false`. |
| `pool` | no | Library pool ID, used only for the read-count scaled bigWigs. Default: one pool for all samples. |

The pipeline checks the sheet before any job starts. It lists every problem at
once and stops if it finds one of these:

- the same sample name with different `control`, `is_control` or `pool` values
- a `control` that is not in the sheet, or that does not have `is_control = true`
- a control sample that has a `control` itself
- a FASTQ file that does not exist or cannot be read, or one listed twice
- a row without both `fastq_1` and `fastq_2`

To check a sheet without running anything:

```bash
nextflow run jacob-greene/CnT-processing --validate_only --genome hg38 --samplesheet samplesheet.csv
```

## IgG controls

- Controls go through the same steps as targets: alignment, QC, bigWigs,
  featureCounts and fragment sizes. They appear next to the targets in every
  table and plot.
- Controls are not peak-called.
- A target with a `control` gets SEACR `<target bedGraph> <control bedGraph> norm stringent`.
  These are its primary peaks, used for FRiP and the plots. The pipeline also
  runs the numeric mode, `<target bedGraph> 0.01 non stringent`, and saves it in
  `seacr/top0.01/` for comparison.
- A target without a `control` gets only the numeric mode (top 1% of signal).
- One control can serve many targets. Each target uses only its own control.
- `control_frip` = control fragments in the target's peaks / all control
  fragments. It is reported next to the target's FRiP.
- The peak plots of each target show its control as an extra track.
- Duplicates are kept for controls too. With `dedup_controls = true`, SEACR uses
  deduplicated control fragments. The Henikoff lab tutorial suggests this,
  because IgG libraries have a high duplication rate. Every other step still
  uses all control fragments.
- If no sample has `is_control = true`, the pipeline runs in numeric mode only
  and says so in the log.

## Supported genomes

Set `genome` to `hg38`, `hg19`, `mm10` or `dm6`. Each genome has pinned
references in `conf/genomes.config`:

| Genome | Annotation (GTF) | cCREs | Blacklist |
|---|---|---|---|
| hg38 | GENCODE 50 | ENCODE v4, ENCFF420VPZ | ENCODE blacklist v2 |
| hg19 | GENCODE 19 | ENCODE v1, ENCFF788SJC (ENCODE has no newer hg19 release) | ENCODE blacklist v2 |
| mm10 | GENCODE M25 | ENCODE v4, ENCFF167FJQ | ENCODE blacklist v2 |
| dm6 | Ensembl 116 (BDGP6.54), renamed to UCSC chromosome names | none: cCRE counting is skipped | ENCODE blacklist v2 |

**Canonical chromosomes.** The UCSC `chrom.sizes` file is used without every
contig whose name has `_` in it (`_random`, `chrUn_*`, `_alt`, `_fix`, `_hap`).
`chrM` is kept (`keep_chrM = true`), in the fragments and in `genome_size` alike.

**genome_size** is computed at run time as the sum of the canonical chromosome
sizes. It is never typed in. With `blacklist_filter = true`, the merged blacklist
bp are subtracted.

| Genome | keep_chrM | genome_size | Blacklist bp | genome_size with blacklist filter |
|---|---|---|---|---|
| hg38 | true | 3,088,286,401 | 227,162,400 | 2,861,124,001 |
| hg38 | false | 3,088,269,832 | 227,162,400 | 2,861,107,432 |
| hg19 | true | 3,095,693,983 | 274,970,000 | 2,820,723,983 |
| hg19 | false | 3,095,677,412 | 274,970,000 | 2,820,707,412 |
| mm10 | true | 2,725,537,669 | 238,977,200 | 2,486,560,469 |
| mm10 | false | 2,725,521,370 | 238,977,200 | 2,486,544,170 |
| dm6 | true | 137,567,484 | 3,740,300 | 133,827,184 |
| dm6 | false | 137,547,960 | 3,740,300 | 133,807,660 |

**Bowtie2 index.** On the Fred Hutch cluster, the default is the shared
iGenomes (UCSC) index of each genome. Elsewhere, set `bowtie2_index` to the
prefix of your own index (the path without `.1.bt2`). To build one from a
FASTA, set `bowtie2_index = false` and `fasta` to the FASTA file(s). The
chromosome names must be UCSC names (`chr1`, ...).

**E. coli index.** Built at each run from `assets/U00096.3.fa.gz` (E. coli
K-12 MG1655, shipped with the pipeline). To skip the build, set `ecoli_index`
to a prebuilt Bowtie2 index prefix.

## Normalization

`C_i` is the total per-base coverage of sample i: the sum of its fragment
lengths. `genome_size` is defined [above](#supported-genomes).

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

## Blacklist filtering (optional)

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

## Parameters

Parameters go in a JSON file (`params.json`), or on the command line as
`--name value`. `examples/params.json` is a complete example:

```json
{
    "samplesheet": "examples/samplesheet.csv",
    "genome": "hg38",
    "outdir": "results",
    "readcount_scaled_bw": true,
    "blacklist_filter": false
}
```

| Parameter | Default | Meaning |
|---|---|---|
| `samplesheet` | required | Sample sheet CSV (see [Sample sheet](#sample-sheet)) |
| `genome` | required | `hg38`, `hg19`, `mm10` or `dm6` |
| `outdir` | `results` | Output folder |
| `min_mapq` | `10` | Minimum MAPQ of a kept read |
| `promoter_window` | `5000` | Promoters = protein-coding TSS +/- this many bp |
| `plot_window` | `5000` | +/- bp around TSSs and peak summits in the deepTools plots |
| `optical_pixel_distance` | `2500` | Picard `OPTICAL_DUPLICATE_PIXEL_DISTANCE`. Use 2500 for patterned flow cells (NovaSeq, HiSeq X/3000/4000) and 100 for unpatterned flow cells (HiSeq 2500, NextSeq, MiSeq) |
| `readcount_scaled_bw` | `false` | Also write read-count scaled bigWigs |
| `dedup_controls` | `false` | SEACR uses deduplicated control fragments |
| `blacklist_filter` | `false` | Remove blacklisted fragments (see above) |
| `keep_chrM` | `true` | Keep `chrM` in the fragments, bigWigs and `genome_size` |
| `bowtie2_index` | per genome | Bowtie2 index prefix; `false` = build it from `fasta` |
| `chrom_sizes`, `gtf`, `ccre_bed`, `blacklist`, `fasta`, `gtf_chrom_map`, `mito_contig` | per genome | Replace one reference of the selected genome |
| `ecoli_fasta`, `ecoli_index` | shipped FASTA | E. coli reference for the carry-over QC |
| `validate_only` | `false` | Check the sample sheet and settings, then stop |
| `max_cpus`, `max_memory`, `max_time` | `16`, `120.GB`, `24.h` | Upper limits for any one job |

**Optical duplicates** need Illumina read names (with the flow-cell
coordinates). Reads renamed by SRA (for example `SRR12246717.1`) have none, so
the optical duplicate columns are `NA` and the log has a warning.

## Launching with `nextflow run`

You need [Nextflow](https://www.nextflow.io/) 24.04 or later, and Docker,
Singularity / Apptainer, or (at Fred Hutch) the Lmod modules.

```bash
nextflow run jacob-greene/CnT-processing -profile <profile> -params-file params.json
```

| Profile | Use |
|---|---|
| `docker` | Pinned biocontainers with Docker |
| `singularity` | The same containers with Singularity / Apptainer |
| `slurm` | Fred Hutch cluster: one Slurm job per task, Singularity containers. Set `--scratch_dir` |
| `modules` | Fred Hutch Lmod modules instead of containers. Put it last, for example `-profile slurm,modules` |
| `test` | The small chr22 + chrM test data set. Combine it with one of the above |

On the Fred Hutch cluster, run the Nextflow head job inside a Slurm job or a
`grabnode` session, for example:

```bash
ml Nextflow
nextflow run jacob-greene/CnT-processing -profile slurm -params-file params.json \
    --scratch_dir /fh/scratch/delete90/<lab>/<user>/nextflow/
```

## Quickstart with the BASH Workbench

The [BASH Workbench](https://github.com/FredHutch/bash-workbench/wiki) is a
menu-driven tool to set up and launch workflows. It can be installed with
`pip3 install bash-workbench`, and it is launched with the command `wb`.

> Users of the Fred Hutch computing cluster can launch the workbench directly
> with `wb`, without any installation.

### Setup

To set up this workflow in the BASH Workbench, do this once:

- Launch the BASH Workbench (`wb`);
- Select `Manage Repositories`;
- Select `Download New Repository`;
- Enter `jacob-greene/CnT-processing` and confirm.

After the setup, exit the workbench with Control+C.

### Launching the workflow

Each run writes its results into the folder that the workbench was launched
from:

- Make a new folder for the output files and go into it (`mkdir my_run && cd my_run`);
- Launch the BASH Workbench (`wb`);
- Select `Run Tool`;
- Select `jacob-greene_CnT-processing`;
- Select `CnT-processing`;
- Enter [the parameters](#parameters). `samplesheet` and `genome` are required;
  the others have defaults;
- Select `Review and Run`;
- Select `jacob-greene_CnT-processing`;
- Select `slurm` (on the Fred Hutch cluster) or `docker` (on a computer with Docker);
- Enter the launcher parameters. For `slurm`, set `scratch_dir` to a folder on the
  scratch file system for temporary files, for example
  `/fh/scratch/delete90/<lab>/<user>/nextflow/`;
- Select `Review and Run`;
- Select `Run Now`.

The workbench keeps a record of the parameters and every log of the run in the
`._wb/` folder of the output folder.

## Resuming and cleaning up

**Resume.** If a run stops (for example, a job ran out of time), start it
again with `-resume` added. Finished steps are not run again. The workbench
`slurm` and `docker` launchers always add `-resume`, so in the workbench, run
the tool again from the same folder.

**Clean up.** Nextflow keeps every intermediate file in its `work/` folder. With
the `slurm` profile or launcher, this folder is `<scratch_dir>/work/`. When the
run is finished and you have checked the results, delete it:

```bash
rm -rf /fh/scratch/delete90/<lab>/<user>/nextflow/work/
```

After this, `-resume` starts from the beginning. Files on `/fh/scratch/delete90`
are deleted automatically 90 days after their last change.

## Output files

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

## Testing

```bash
nextflow run jacob-greene/CnT-processing -profile test,docker --outdir results
python3 tests/check_outputs.py results --test-data    # needs numpy and pyBigWig
```

The test data are three K562 samples from Kaya-Okur et al. 2020 (GEO GSE145187):
H3K27me3, H3K4me3 and IgG, restricted to read pairs on chr22 and chrM, at three
different depths. IgG is split into two lanes. Every test input is in
`tests/data/`, so the test run downloads nothing. `tests/make_test_data.sh`
rebuilds the test data and also writes the larger smoke-test data (1M read pairs
per sample) to `test_data_work/smoke/`. `examples/params.json` and
`examples/samplesheet.csv` run the smoke test against full hg38.

`tests/check_outputs.py` checks known answers, for example:

- every coverage-normalized bigWig averages 1 over the genome, and on 80 random
  windows each base equals the fragment pileup times its scale factor. A
  negative control shuffles bigWig values between equal-length intervals (same
  total); the per-base check must then fail
- the read-count scaled bigWigs keep the coverage ratios between samples
- every featureCounts value equals the number of fragments that overlap that
  feature, counted again from the fragment and SAF files
- FRiP, control FRiP, peak counts and widths, recomputed from the fragment and
  peak files
- the 5 kb bin means, rebuilt from the fragment files, equal the pipeline's bin
  file and the bigWigs; the correlation matrix is recomputed from these bins

GitHub Actions runs the test profile twice on every push, with
`blacklist_filter = false` and `true`, then the checks.

## Authors

Jacob Greene

## License

MIT
