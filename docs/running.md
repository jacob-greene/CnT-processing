# Running the pipeline

- [Parameters](#parameters)
- [Launching with `nextflow run`](#launching-with-nextflow-run)
- [Quickstart with the BASH Workbench](#quickstart-with-the-bash-workbench)
- [Resuming and cleaning up](#resuming-and-cleaning-up)

Back to the [README](../README.md).

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
| `samplesheet` | required | Sample sheet CSV (see [Sample sheet](inputs.md#sample-sheet)) |
| `genome` | required | `hg38`, `hg19`, `mm10` or `dm6` |
| `outdir` | `results` | Output folder |
| `min_mapq` | `10` | Minimum MAPQ of a kept read |
| `promoter_window` | `5000` | Promoters = protein-coding TSS +/- this many bp |
| `plot_window` | `5000` | +/- bp around TSSs and peak summits in the deepTools plots |
| `optical_pixel_distance` | `2500` | Picard `OPTICAL_DUPLICATE_PIXEL_DISTANCE`. Use 2500 for patterned flow cells (NovaSeq, HiSeq X/3000/4000) and 100 for unpatterned flow cells (HiSeq 2500, NextSeq, MiSeq) |
| `readcount_scaled_bw` | `false` | Also write read-count scaled bigWigs |
| `dedup_controls` | `false` | SEACR uses deduplicated control fragments |
| `blacklist_filter` | `false` | Remove blacklisted fragments (see [Blacklist filtering](steps.md#blacklist-filtering-optional)) |
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
