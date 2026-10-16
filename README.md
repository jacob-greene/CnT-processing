# CnT-processing

A small, readable Nextflow (DSL2) pipeline for paired-end CUT&Tag data.
It takes FASTQ files and produces:

1. Coverage-normalized bigWigs
2. Optional read-count scaled bigWigs
3. featureCounts matrices over protein-coding genes, promoters and ENCODE cCREs
4. SEACR peaks (IgG-controlled when a control is given), FRiP, peak stats and deepTools plots
5. QC tables: alignment stats, % chrM, fragment sizes, duplication, E. coli carry-over, sample correlation
6. A MultiQC report

> **Status: under construction (Phase 1 of 4).** The sample sheet validator,
> the genome configuration and the test-data builder are in place. The
> processing steps are not implemented yet. The full user guide, including the
> BASH Workbench quickstart, is written in Phase 4.

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
| `bin/` | Scripts used by the processes |
| `assets/` | Small static files (e.g. the dm6 chromosome name map) |
| `tests/` | Sample sheet validator tests and `make_test_data.sh` |
| `._wb/` | BASH Workbench tool and launcher definitions |

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
