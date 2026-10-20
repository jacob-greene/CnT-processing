# Inputs

- [Input data](#input-data)
- [Sample sheet](#sample-sheet)
- [IgG controls](#igg-controls)
- [Supported genomes](#supported-genomes)

Back to the [README](../README.md).

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
