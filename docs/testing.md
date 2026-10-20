# Testing

Back to the [README](../README.md).

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
