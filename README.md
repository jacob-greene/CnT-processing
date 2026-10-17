# CnT-processing

A small, readable Nextflow (DSL2) pipeline for paired-end CUT&Tag data.
It takes FASTQ files and produces:

1. Coverage-normalized bigWigs
2. Optional read-count scaled bigWigs
3. featureCounts matrices over protein-coding genes, promoters and ENCODE cCREs
4. SEACR peaks (IgG-controlled when a control is given), FRiP, peak stats and deepTools plots
5. QC tables: alignment stats, % chrM, fragment sizes, duplication, E. coli carry-over, sample correlation
6. A MultiQC report

> **Status: under construction (Phase 2 of 4).** Steps 1-10 below are
> implemented and tested: FASTQ to filtered BAMs, QC tables, fragment files and
> both kinds of bigWig. Steps 11-16 (featureCounts, SEACR peaks, FRiP,
> deepTools plots, correlation, MultiQC) come in Phase 3. The full user guide,
> including the BASH Workbench quickstart, is written in Phase 4.

## Repository layout

| Path | Contents |
|---|---|
| `main.nf` | Workflow entry point |
| `modules/` | Processes and helper functions (sample sheet parser, genome preparation, ...) |
| `nextflow.config` | Parameters, containers and profiles (`docker`, `singularity`, `slurm`, `test`) |
| `conf/genomes.config` | Pinned reference URLs and paths for `hg38`, `hg19`, `mm10`, `dm6` |
| `conf/base.config` | Default resources and out-of-memory retry |
| `conf/slurm.config` | Fred Hutch gizmo cluster settings |
| `conf/test.config` | Small chr22 + chrM test data set |
| `conf/modules.config` | Fred Hutch Lmod modules instead of containers (`-profile modules`) |
| `bin/` | Scripts used by the processes |
| `assets/` | Small static files (e.g. the dm6 chromosome name map) |
| `tests/` | Validator tests, `make_test_data.sh`, the test data and `check_outputs.py` |
| `.github/workflows/ci.yml` | GitHub Actions: the test profile with Docker, blacklist filter off and on |
| `._wb/` | BASH Workbench tool and launcher definitions |

## Workflow steps (implemented so far)

| Step | What happens | Tool |
|---|---|---|
| 1 | Join the lanes of each sample, then run FastQC | cat, FastQC |
| 2 | Adapter and quality trimming (pairs shorter than 20 bp are dropped) | Trim Galore (cutadapt) |
| 3 | Align to the genome: `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700` | Bowtie2 |
| 4 | E. coli carry-over QC: the same reads against E. coli K-12 MG1655, plus `--no-overlap --no-dovetail`. QC only, no spike-in normalization | Bowtie2 |
| 5 | Keep properly paired, primary pairs with MAPQ >= `min_mapq` (`-f 2 -F 0x904 -q`); drop a mate whose partner failed the filter; sort and index | samtools |
| 6 | Mark duplicates, but keep them (`REMOVE_DUPLICATES=false`) | Picard MarkDuplicates |
| 7 | `alignment_stats.tsv`: one row per sample | Python |
| 8 | Fragment size counts (`fragment_sizes/<sample>.tsv`) and one combined plot | samtools, matplotlib |
| 9 | Fragment BED (`bedtools bamtobed -bedpe`; same chromosome, < 1000 bp, canonical chromosomes); optional blacklist filter | bedtools |
| 10 | Raw-count bedGraphs (for SEACR), `scale_factors.tsv`, coverage-normalized and (optional) read-count scaled bigWigs | bedtools, bedGraphToBigWig |

## Normalization

`C_i` is the total per-base coverage of sample i, which is the sum of its
fragment lengths. `genome_size` is the sum of the canonical chromosome sizes.
With `blacklist_filter = true`, the merged blacklist bp are subtracted.

**Coverage-normalized bigWigs** (always made, `bigwig/coverage_norm/`):

```
s_cov_i = genome_size / C_i
```

Each base holds the fraction of the sample's total coverage at that base,
times `genome_size`. If coverage were uniform, every base would be 1. The
genome-length-weighted mean of every bigWig is 1.

**Read-count scaled bigWigs** (`readcount_scaled_bw = true`, `bigwig/readcount_scaled/`):

```
s_rc_i  = genome_size / mean_j(C_j)        j = every sample in the same pool
bw_rc_i = bw_cov_i x C_i / mean_j(C_j)
```

The assumption: every library in a pool started from the same number of cells
and was pooled by equal volume. Then differences in fragment yield reflect
epitope abundance, and one shared factor per pool keeps them. The average
sample of a pool averages about 1 per bp. **These bigWigs are only meaningful
if the assumption holds.** Controls are included in their pool.

`scale_factors.tsv` lists, per sample: pool, fragments, `C_i`, `s_cov`,
`s_rc`, `C_i / mean_j(C_j)`, the blacklist filter state, the fragments removed
by the blacklist and the `genome_size` used.

## Profiles

| Profile | Use |
|---|---|
| `docker` | Pinned biocontainers with Docker |
| `singularity` | The same containers with Singularity / Apptainer |
| `slurm` | Fred Hutch gizmo: Slurm jobs + Singularity; set `--scratch_dir` |
| `modules` | Fred Hutch Lmod modules instead of containers; put it last, e.g. `-profile slurm,modules` |
| `test` | The small chr22 + chrM test data; combine with one of the above |

## Test run

```bash
nextflow run jacob-greene/CnT-processing -profile test,docker --outdir results
python3 tests/check_outputs.py results --test-data    # needs pyBigWig
```

`tests/check_outputs.py` checks the known answers: every coverage-normalized
bigWig averages 1 over the genome, the read-count scaled bigWigs keep the
coverage ratios between samples, the scale factors match a recomputation from
the fragment files, and the alignment table columns reconcile.

Without `--ecoli_index`, each run builds the E. coli Bowtie2 index from
`assets/U00096.3.fa.gz` (ENA U00096.3, shipped with the pipeline). To skip the
build, pass a prebuilt Bowtie2 index with `--ecoli_index`.

## Check a sample sheet

```bash
nextflow run jacob-greene/CnT-processing --validate_only \
    --genome hg38 --samplesheet samplesheet.csv
```

## Sample sheet

```
sample,fastq_1,fastq_2,control,is_control,pool
H3K27me3_rep1,/path/R1.fastq.gz,/path/R2.fastq.gz,IgG_rep1,false,pool1
H3K4me3_rep1,/path/R1.fastq.gz,/path/R2.fastq.gz,IgG_rep1,false,pool1
IgG_rep1,/path/R1.fastq.gz,/path/R2.fastq.gz,,true,pool1
```

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample name. Rows with the same name are lanes and are concatenated. |
| `fastq_1`, `fastq_2` | yes | R1 and R2 FASTQ. A relative path is resolved against the sample sheet's folder. |
| `control` | no | Name of an IgG / negative-control sample in the same sheet. |
| `is_control` | no | `true` for IgG / negative controls. Default `false`. |
| `pool` | no | Library pool ID for read-count scaling. Default: one shared pool. |

The pipeline stops before any job starts if the sheet has a problem, and it
lists every problem at once.

## Authors

Jacob Greene

## License

MIT
