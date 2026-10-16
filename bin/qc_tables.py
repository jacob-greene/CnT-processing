#!/usr/bin/env python3
"""Build alignment_stats.tsv and duplication.tsv (one row per sample).

All per-sample files are read from the current folder, by name:
    <id>.read_pairs.tsv          total_pairs (from CAT_FASTQ)
    <id>.bowtie2.log             Bowtie2 log, target genome
    <id>.ecoli.bowtie2.log       Bowtie2 log, E. coli
    <id>.filter.tsv              filtered_pairs, mito_pairs, illumina_read_names
    <id>.markdup_metrics.txt     Picard MarkDuplicates metrics
    <id>.fragments.tsv           fragments (final, after any blacklist filter)

Usage: qc_tables.py id1 id2 ...
"""
import re
import sys


def read_kv(path):
    with open(path) as fh:
        return dict(line.rstrip('\n').split('\t', 1) for line in fh if '\t' in line)


def bowtie2_pairs(path):
    """Return (input pairs, concordant pairs) from a paired-end Bowtie2 log."""
    text = open(path).read()
    total = int(re.search(r'^(\d+) reads; of these:', text, re.M).group(1))
    once = int(re.search(r'(\d+) \([\d.]+%\) aligned concordantly exactly 1 time', text).group(1))
    multi = int(re.search(r'(\d+) \([\d.]+%\) aligned concordantly >1 times', text).group(1))
    return total, once + multi


def picard_metrics(path):
    lines = open(path).read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith('## METRICS CLASS'):
            header = lines[i + 1].split('\t')
            values = lines[i + 2].split('\t')
            return dict(zip(header, values))
    sys.exit(f'ERROR: no metrics in {path}')


def pct(num, den):
    return f'{100.0 * num / den:.2f}' if den else 'NA'


def main(samples):
    align_cols = ['sample', 'total_pairs', 'pairs_after_trimming', 'aligned_pairs', 'alignment_rate_pct',
                  'filtered_pairs', 'pct_chrM', 'duplication_pct', 'optical_duplicate_pct',
                  'ecoli_pairs', 'ecoli_pct', 'final_fragments']
    dup_cols = ['sample', 'total_pairs', 'duplicate_pairs', 'PERCENT_DUPLICATION',
                'optical_duplicate_pairs', 'optical_duplicate_fraction', 'ESTIMATED_LIBRARY_SIZE']
    align_rows, dup_rows, warnings = [], [], []

    for s in samples:
        total = int(read_kv(f'{s}.read_pairs.tsv')['total_pairs'])
        trimmed, aligned = bowtie2_pairs(f'{s}.bowtie2.log')
        ecoli_in, ecoli = bowtie2_pairs(f'{s}.ecoli.bowtie2.log')
        filt = read_kv(f'{s}.filter.tsv')
        filtered = int(filt['filtered_pairs'])
        mito = int(filt['mito_pairs'])
        illumina = filt['illumina_read_names'] == 'true'
        m = picard_metrics(f'{s}.markdup_metrics.txt')
        frags = read_kv(f'{s}.fragments.tsv')['fragments']

        pairs = int(m['READ_PAIRS_EXAMINED'])
        dup_pairs = int(m['READ_PAIR_DUPLICATES'])
        dup_pct = 100.0 * float(m['PERCENT_DUPLICATION'] or 0)
        if illumina:
            optical = int(m['READ_PAIR_OPTICAL_DUPLICATES'])
            optical_frac = f'{optical / pairs:.6f}' if pairs else 'NA'
            optical_pct = pct(optical, pairs)
        else:
            optical, optical_frac, optical_pct = 'NA', 'NA', 'NA'
            warnings.append(f'{s}: read names have no Illumina flow-cell coordinates '
                            f'(e.g. "{filt["first_read_name"]}"); optical duplicates reported as NA')

        align_rows.append([s, total, trimmed, aligned, pct(aligned, trimmed), filtered,
                           pct(mito, filtered), f'{dup_pct:.2f}', optical_pct,
                           ecoli, pct(ecoli, ecoli_in), frags])
        dup_rows.append([s, pairs, dup_pairs, m['PERCENT_DUPLICATION'], optical, optical_frac,
                         m['ESTIMATED_LIBRARY_SIZE'] or 'NA'])

        # The columns must reconcile.
        if not (filtered <= aligned <= trimmed <= total):
            sys.exit(f'ERROR: {s}: pair counts do not reconcile: filtered {filtered}, '
                     f'aligned {aligned}, after trimming {trimmed}, total {total}')
        if pairs != filtered:
            sys.exit(f'ERROR: {s}: Picard examined {pairs} pairs, but {filtered} pairs passed the filter')

    for name, cols, rows in [('alignment_stats.tsv', align_cols, align_rows),
                             ('duplication.tsv', dup_cols, dup_rows)]:
        with open(name, 'w') as fh:
            fh.write('\t'.join(cols) + '\n')
            for r in rows:
                fh.write('\t'.join(str(x) for x in r) + '\n')
    for w in warnings:
        print(f'WARNING: {w}', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1:])
