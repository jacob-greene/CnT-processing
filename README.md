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

- [Workflow](#workflow)
- [Quick start](#quick-start)
- [Sample sheet](#sample-sheet)
- [Key parameters](#key-parameters)
- [Main outputs](#main-outputs)
- [More documentation](#more-documentation)

## Workflow

```mermaid
flowchart TD
    fq[/"FASTQ lanes + sample sheet"/]
    ref[/"Genome: GENOME_SIZE, PREPARE_GTF, BOWTIE2_BUILD"/]

    fq --> CAT_FASTQ["CAT_FASTQ (1)"]
    CAT_FASTQ --> FASTQC["FASTQC (1)"]
    CAT_FASTQ --> TRIM_GALORE["TRIM_GALORE (2)"]
    TRIM_GALORE --> BOWTIE2_ALIGN["BOWTIE2_ALIGN (3)"]
    TRIM_GALORE --> BOWTIE2_ECOLI["BOWTIE2_ECOLI (4)<br/>QC only"]
    ref --> BOWTIE2_ALIGN
    ref --> BOWTIE2_ECOLI
    BOWTIE2_ALIGN --> FILTER_BAM["FILTER_BAM (5)"]
    FILTER_BAM --> MARKDUPLICATES["MARKDUPLICATES (6)<br/>duplicates kept"]
    MARKDUPLICATES --> NAMESORT
    NAMESORT --> FRAGMENTS["FRAGMENTS (9)<br/>fragment BED + raw bedGraph<br/>blacklist filter if blacklist_filter"]
    NAMESORT -. "dedup_controls" .-> FRAGMENTS_DEDUP["FRAGMENTS_DEDUP<br/>deduplicated IgG"]

    FRAGMENTS --> SCALE_FACTORS["SCALE_FACTORS (10)"]
    SCALE_FACTORS --> BIGWIG_COV["BIGWIG_COV (10)<br/>coverage-normalized"]
    SCALE_FACTORS -. "readcount_scaled_bw" .-> BIGWIG_RC["BIGWIG_RC (10)<br/>read-count scaled"]

    ref --> ANNOTATIONS["ANNOTATIONS<br/>genes, promoters, cCREs, TSSs"]
    FRAGMENTS --> FEATURECOUNTS["FEATURECOUNTS (11)"]
    ANNOTATIONS --> FEATURECOUNTS

    FRAGMENTS --> SEACR{"SEACR (12)<br/>targets only"}
    FRAGMENTS_DEDUP -. "control bedGraph" .-> SEACR
    SEACR -- "control set" --> ctl["IgG-controlled<br/>norm stringent, + top 1% set"]
    SEACR -- "no control" --> noctl["top 1%<br/>0.01 non stringent"]
    ctl --> PEAK_STATS["PEAK_STATS (13)<br/>FRiP, control FRiP"]
    noctl --> PEAK_STATS
    ctl --> DEEPTOOLS_PEAKS["DEEPTOOLS_PEAKS (14)"]
    noctl --> DEEPTOOLS_PEAKS
    BIGWIG_COV --> DEEPTOOLS_PEAKS
    BIGWIG_COV --> DEEPTOOLS_TSS["DEEPTOOLS_TSS (14)"]
    BIGWIG_RC -.-> DEEPTOOLS_TSS
    ANNOTATIONS --> DEEPTOOLS_TSS
    BIGWIG_COV --> CORRELATION["CORRELATION (15)<br/>5 kb bins"]

    subgraph qc["QC"]
        QC_TABLES["QC_TABLES (7)<br/>alignment stats, duplication"]
        FRAGMENT_SIZE_PLOT["FRAGMENT_SIZE_PLOT (8)"]
        MULTIQC["MULTIQC (16)"]
        QC_TABLES --> MULTIQC
        FRAGMENT_SIZE_PLOT --> MULTIQC
    end
    BOWTIE2_ECOLI --> QC_TABLES
    FILTER_BAM --> FRAGMENT_SIZE_PLOT
    FASTQC --> MULTIQC

    BIGWIG_COV --> o1[("bigWigs")]
    BIGWIG_RC -.-> o1
    FEATURECOUNTS --> o2[("featureCounts matrices")]
    ctl --> o3[("SEACR peaks")]
    noctl --> o3
    PEAK_STATS --> o4[("peaks_frip.tsv")]
    DEEPTOOLS_PEAKS --> o5[("deepTools plots")]
    DEEPTOOLS_TSS --> o5
    CORRELATION --> o5
    MULTIQC --> o6[("MultiQC report")]

    classDef io fill:#e0e0e0,stroke:#616161,color:#000
    classDef prep fill:#d7b899,stroke:#795548,color:#000
    classDef align fill:#9ecae1,stroke:#2171b5,color:#000
    classDef frag fill:#cbb4e8,stroke:#6a3d9a,color:#000
    classDef cov fill:#a6dba0,stroke:#1b7837,color:#000
    classDef counts fill:#fee391,stroke:#b8860b,color:#000
    classDef peaks fill:#fdbe85,stroke:#d94801,color:#000
    classDef qcs fill:#f4a6a6,stroke:#c62828,color:#000

    class fq,ref,o1,o2,o3,o4,o5,o6 io
    class CAT_FASTQ,FASTQC,TRIM_GALORE prep
    class BOWTIE2_ALIGN,BOWTIE2_ECOLI,FILTER_BAM,MARKDUPLICATES,NAMESORT align
    class FRAGMENTS,FRAGMENTS_DEDUP frag
    class SCALE_FACTORS,BIGWIG_COV,BIGWIG_RC cov
    class ANNOTATIONS,FEATURECOUNTS counts
    class SEACR,ctl,noctl,PEAK_STATS,DEEPTOOLS_PEAKS,DEEPTOOLS_TSS,CORRELATION peaks
    class QC_TABLES,FRAGMENT_SIZE_PLOT,MULTIQC qcs
    style qc fill:transparent,stroke:#c62828,stroke-dasharray:4 4
```

Each color is one stage. The Details column links to the full description.

| Color | Stage | Steps | Details |
|---|---|---|---|
| ⬜ grey | Inputs and outputs | - | [Inputs](docs/inputs.md), [Output files](docs/outputs.md) |
| 🟫 brown | Read preparation | 1-2 | [Read preparation](docs/steps.md#read-preparation) |
| 🟦 blue | Alignment and filtering | 3-6 | [Alignment and filtering](docs/steps.md#alignment-and-filtering) |
| 🟪 purple | Fragments | 9 | [Fragments](docs/steps.md#fragments) |
| 🟩 green | Coverage and bigWigs | 10 | [Coverage and bigWigs](docs/steps.md#coverage-and-bigwigs) |
| 🟨 yellow | Feature counts | 11 | [Feature counts](docs/steps.md#feature-counts) |
| 🟧 orange | Peaks and signal plots | 12-15 | [Peaks and signal plots](docs/steps.md#peaks-and-signal-plots) |
| 🟥 red | QC reports | 7, 8, 16 | [QC reports](docs/steps.md#qc-reports) |

Solid arrows always run. Dotted arrows run only when the named parameter is
true. Numbers in brackets are step numbers. The [step table](docs/steps.md#all-steps)
lists all 16 steps. Some minor inputs are not drawn; see
[what the diagram leaves out](docs/steps.md#what-the-diagram-leaves-out).

## Quick start

You need [Nextflow](https://www.nextflow.io/) 24.04 or later, and Docker,
Singularity / Apptainer, or (at Fred Hutch) the Lmod modules.

```bash
nextflow run jacob-greene/CnT-processing -profile <profile> -params-file params.json
```

Use the `docker`, `singularity` or (Fred Hutch) `slurm` profile. To try the
pipeline on the small test data set:

```bash
nextflow run jacob-greene/CnT-processing -profile test,docker --outdir results
```

More: [all profiles, the Fred Hutch cluster, the BASH Workbench, resume and
clean-up](docs/running.md), and [the tests](docs/testing.md).

## Sample sheet

A CSV file with one row per FASTQ pair. Rows with the same `sample` are lanes
and are joined.

```
sample,fastq_1,fastq_2,control,is_control,pool
H3K27me3_rep1,/path/to/H3K27me3_rep1_R1.fastq.gz,/path/to/H3K27me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1
IgG_rep1,/path/to/IgG_rep1_R1.fastq.gz,/path/to/IgG_rep1_R2.fastq.gz,,true,pool1
```

`sample`, `fastq_1` and `fastq_2` are required. `control` names the IgG sample
of a target. The pipeline checks the sheet before any job starts. More:
[sample sheet columns and checks](docs/inputs.md#sample-sheet),
[IgG controls](docs/inputs.md#igg-controls) and
[supported genomes](docs/inputs.md#supported-genomes) (`hg38`, `hg19`, `mm10`,
`dm6`).

## Key parameters

Parameters go in a JSON file (`params.json`) or on the command line as
`--name value`. `examples/params.json` is a complete example.

| Parameter | Default | Meaning |
|---|---|---|
| `samplesheet` | required | Sample sheet CSV |
| `genome` | required | `hg38`, `hg19`, `mm10` or `dm6` |
| `outdir` | `results` | Output folder |
| `readcount_scaled_bw` | `false` | Also write [read-count scaled bigWigs](docs/steps.md#coverage-and-bigwigs) |
| `dedup_controls` | `false` | SEACR uses deduplicated control fragments |
| `blacklist_filter` | `false` | Remove [blacklisted fragments](docs/steps.md#blacklist-filtering-optional) |

All parameters: [Parameters](docs/running.md#parameters).

## Main outputs

| Output | What it is |
|---|---|
| `bigwig/coverage_norm/` | Coverage-normalized bigWigs ([how](docs/steps.md#coverage-and-bigwigs)) |
| `featurecounts/` | Fragment counts over genes, promoters and cCREs |
| `seacr/` | SEACR peaks of each target |
| `peaks_frip.tsv` | Peak count, widths, FRiP and control FRiP per target |
| `alignment_stats.tsv` | One QC row per sample |
| `multiqc/multiqc_report.html` | The MultiQC report |

All outputs: [Output files](docs/outputs.md).

## More documentation

| Page | Contents |
|---|---|
| [Workflow steps](docs/steps.md) | Every step, its tool and settings, normalization, blacklist filter, featureCounts and correlation notes |
| [Inputs](docs/inputs.md) | Input data, sample sheet, IgG controls, supported genomes |
| [Running](docs/running.md) | All parameters, profiles, Fred Hutch cluster, BASH Workbench, resume and clean-up |
| [Output files](docs/outputs.md) | The output folder and the `alignment_stats.tsv` columns |
| [Testing](docs/testing.md) | Test data and the output checks |

## Authors

Jacob Greene

## License

MIT
