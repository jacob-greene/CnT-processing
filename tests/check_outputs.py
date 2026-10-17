#!/usr/bin/env python3
"""Known-answer checks on a finished CnT-processing run (steps 1-10).

    python3 tests/check_outputs.py <outdir> [--test-data]

Needs Python 3 with pyBigWig (e.g. the deepTools container or module).
Every check prints PASS or FAIL; the exit status is 1 if any check fails.

Checks:
  * coverage-normalized bigWigs: genome-length-weighted mean = 1
    (total signal / genome_size; with blacklist_filter = true, blacklisted bp
    carry no fragments, so this is the mean over the non-blacklisted genome)
  * read-count scaled bigWigs: mean_i = C_i / mean_j(C_j), pool mean of these = 1,
    and mean_i / mean_k = C_i / C_k for every pair of samples in a pool
  * per-base signal: on a fixed random sample of windows, every bigWig base =
    (fragment pileup rebuilt from the fragment BED) x its scale factor. This
    catches errors that keep the total, e.g. values moved between intervals.
    Negative control (first sample): the same check must fail on a copy of
    each bigWig with values shuffled among equal-length intervals
  * s_cov and s_rc recomputed from the fragment BED files, independent of the pipeline
  * alignment_stats.tsv: every column filled, % chrM > 0 (test data),
    filtered <= aligned <= after trimming <= total
  * duplicates kept: fragments (+ blacklist-removed) = pairs examined by MarkDuplicates
    (--test-data only: every test pair is on chr22/chrM with a fragment < 1000 bp)
  * optical duplicates are NA for the SRA-named test reads (--test-data)
  * blacklist_filter = true: 0 fragments overlap the blacklist, and genome_size =
    canonical bp - merged blacklist bp
"""
import argparse
import bisect
import csv
import gzip
import os
import sys
import tempfile

import numpy as np
import pyBigWig

from mutate_bigwig import shuffle

TOL = 1e-3
WINDOW = 5000           # bp per per-base check window
N_WINDOWS = (60, 20)    # windows centred on random fragments, windows anywhere
failures = 0


def check(ok, what, detail=''):
    global failures
    print(f'{"PASS" if ok else "FAIL"}  {what}{"  (" + detail + ")" if detail else ""}')
    if not ok:
        failures += 1


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def read_kv(path):
    with open(path) as fh:
        return dict(line.rstrip('\n').split('\t', 1) for line in fh if '\t' in line)


def total_signal(path):
    bw = pyBigWig.open(path)
    total = 0.0
    for chrom, length in bw.chroms().items():
        for start, end, value in bw.intervals(chrom) or []:
            total += value * (end - start)
    bw.close()
    return total


def fragments(path):
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            c, s, e = line.split('\t')[:3]
            yield c, int(s), int(e)


def check_windows(frags, sizes, seed=1):
    """Fixed random windows: most centred on fragments (signal), some anywhere (mostly zero)."""
    rng = np.random.default_rng(seed)
    wins = []
    for i in rng.choice(len(frags), size=min(N_WINDOWS[0], len(frags)), replace=False):
        c, b, e = frags[i]
        mid = (b + e) // 2
        wins.append((c, max(0, mid - WINDOW // 2), min(sizes[c], mid + WINDOW // 2)))
    chroms = sorted(sizes)
    weights = np.array([sizes[c] for c in chroms], float)
    for c in rng.choice(chroms, size=N_WINDOWS[1], p=weights / weights.sum()):
        b = int(rng.integers(0, max(1, sizes[c] - WINDOW)))
        wins.append((c, b, min(sizes[c], b + WINDOW)))
    return wins


def pileup(frags_by_chrom, starts_by_chrom, c, b, e):
    """Per-base fragment count on [b, e), rebuilt from the fragments (all < 1000 bp)."""
    cov = np.zeros(e - b + 1)
    fr, st = frags_by_chrom.get(c, []), starts_by_chrom.get(c, [])
    for fb, fe in fr[bisect.bisect_left(st, b - 1000):bisect.bisect_left(st, e)]:
        lo, hi = max(fb, b), min(fe, e)
        if lo < hi:
            cov[lo - b] += 1
            cov[hi - b] -= 1
    return np.cumsum(cov)[:-1]


def per_base_check(sample, frags, sizes, bigwigs, negative_control=False):
    """bigwigs: list of (kind, path, scale factor). One PASS/FAIL per bigWig."""
    by_chrom = {}
    for c, b, e in frags:
        by_chrom.setdefault(c, []).append((b, e))
    starts = {c: [b for b, _ in v] for c, v in by_chrom.items()}
    wins = check_windows(frags, sizes)

    def errors(path, factor):
        bw = pyBigWig.open(path)
        bad, nonzero, worst = 0, 0, 0.0
        for c, b, e in wins:
            want = pileup(by_chrom, starts, c, b, e) * factor
            got = np.nan_to_num(np.array(bw.values(c, b, e), dtype=float))
            err = np.abs(got - want) / np.maximum(want, 1e-9)
            err[want == 0] = np.abs(got[want == 0])
            bad += int((err > 1e-4).sum())
            worst = max(worst, float(err.max()))
            nonzero += int((want > 0).sum())
        bw.close()
        return bad, nonzero, worst

    for kind, path, factor in bigwigs:
        bad, nonzero, worst = errors(path, factor)
        check(bad == 0 and nonzero > 0, f'{sample}: {kind} per-base = fragment pileup x factor',
              f'{len(wins)} windows, {nonzero:,} covered bp, {bad} bp off, max rel err {worst:.1e}')
        if negative_control:
            # The same check on a copy with values shuffled among equal-length
            # intervals (total unchanged) must fail, or it cannot catch such errors.
            with tempfile.TemporaryDirectory() as tmp:
                mutated = os.path.join(tmp, 'mutated.bw')
                moved = shuffle(path, mutated)
                bad_m = errors(mutated, factor)[0]
                d_total = abs(total_signal(mutated) - total_signal(path)) / total_signal(path)
            check(moved > 0 and bad_m > 0 and d_total < 1e-6,
                  f'{sample}: negative control - per-base check fails on shuffled {kind} bigWig',
                  f'{moved:,} values moved, {bad_m:,} bp off, total changed by {d_total:.1e}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('--test-data', action='store_true',
                    help='also run the checks that only hold for the chr22 + chrM test data')
    args = ap.parse_args()
    out = args.outdir

    gs = read_kv(f'{out}/genome/genome_size.tsv')
    genome_size = int(gs['genome_size'])
    blacklist_on = gs['blacklist_filter'] == 'true'
    sf = read_tsv(f'{out}/scale_factors.tsv')
    stats = {r['sample']: r for r in read_tsv(f'{out}/alignment_stats.tsv')}
    dup = {r['sample']: r for r in read_tsv(f'{out}/duplication.tsv')}
    print(f'genome_size = {genome_size:,} bp, blacklist_filter = {gs["blacklist_filter"]}\n')

    # ---- Independent recomputation of C_i, s_cov, s_rc from the fragment BEDs ----
    C = {}
    for r in sf:
        s = r['sample']
        C[s] = sum(e - b for _, b, e in fragments(f'{out}/fragments/{s}.fragments.bed.gz'))
        check(C[s] == int(r['C_i']), f'{s}: C_i from fragment BED = scale_factors.tsv', f'{C[s]:,}')
        s_cov = genome_size / C[s]
        check(abs(s_cov / float(r['s_cov']) - 1) < 1e-6, f'{s}: s_cov recomputed', f'{s_cov:.6f}')
    pools = {}
    for r in sf:
        pools.setdefault(r['pool'], []).append(r['sample'])
    for r in sf:
        members = pools[r['pool']]
        mean = sum(C[m] for m in members) / len(members)
        check(abs((genome_size / mean) / float(r['s_rc']) - 1) < 1e-6, f'{r["sample"]}: s_rc recomputed',
              f'{genome_size / mean:.6f}')
        check(r['genome_size'] == str(genome_size), f'{r["sample"]}: genome_size column matches genome_size.tsv')

    # ---- Coverage-normalized bigWigs ----
    print()
    for r in sf:
        s = r['sample']
        m = total_signal(f'{out}/bigwig/coverage_norm/{s}.coverage_norm.bw') / genome_size
        check(abs(m - 1) < TOL, f'{s}: coverage_norm weighted mean = 1', f'{m:.6f}')

    # ---- Per-base signal on a sample of windows ----
    print()
    sizes = {}
    with open(f'{out}/genome/canonical.chrom.sizes') as fh:
        for line in fh:
            c, n = line.split('\t')[:2]
            sizes[c] = int(n)
    rc_dir = f'{out}/bigwig/readcount_scaled'
    for r in sf:
        s = r['sample']
        bigwigs = [('coverage_norm', f'{out}/bigwig/coverage_norm/{s}.coverage_norm.bw', float(r['s_cov']))]
        if os.path.isdir(rc_dir):
            bigwigs.append(('readcount_scaled', f'{rc_dir}/{s}.readcount_scaled.bw', float(r['s_rc'])))
        per_base_check(s, list(fragments(f'{out}/fragments/{s}.fragments.bed.gz')), sizes, bigwigs,
                       negative_control=(r is sf[0]))

    # ---- Read-count scaled bigWigs ----
    if os.path.isdir(rc_dir):
        print()
        means = {}
        for r in sf:
            s = r['sample']
            means[s] = total_signal(f'{rc_dir}/{s}.readcount_scaled.bw') / genome_size
            ratio = float(r['C_i_over_pool_mean'])
            check(abs(means[s] - ratio) < TOL * max(1, ratio), f'{s}: readcount_scaled mean = C_i / mean_j(C_j)',
                  f'{means[s]:.6f} vs {ratio:.6f}')
        for pool, members in pools.items():
            pm = sum(means[m] for m in members) / len(members)
            check(abs(pm - 1) < TOL, f'pool {pool}: mean of readcount_scaled means = 1', f'{pm:.6f}')
            for i, a in enumerate(members):
                for b in members[i + 1:]:
                    got, want = means[a] / means[b], C[a] / C[b]
                    check(abs(got / want - 1) < TOL, f'{a} / {b}: bigWig mean ratio = coverage ratio',
                          f'{got:.6f} vs {want:.6f}')
    else:
        print('\n(no readcount_scaled bigWigs: readcount_scaled_bw = false)')

    # ---- Alignment stats ----
    print()
    for s, r in stats.items():
        empty = [k for k, v in r.items() if v in ('', None) or (v == 'NA' and k != 'optical_duplicate_pct')]
        check(not empty, f'{s}: every alignment_stats column filled', ', '.join(empty))
        tot, trim, al, filt = (int(r[k]) for k in ('total_pairs', 'pairs_after_trimming', 'aligned_pairs',
                                                   'filtered_pairs'))
        check(filt <= al <= trim <= tot, f'{s}: filtered <= aligned <= after trimming <= total',
              f'{filt} <= {al} <= {trim} <= {tot}')
        if args.test_data:
            check(float(r['pct_chrM']) > 0, f'{s}: % chrM > 0', r['pct_chrM'])

    # ---- Duplicates kept ----
    print()
    for r in sf:
        s = r['sample']
        examined = int(dup[s]['total_pairs'])
        removed = 0 if r['fragments_removed_blacklist'] == 'NA' else int(r['fragments_removed_blacklist'])
        kept = int(r['fragments']) + removed
        if args.test_data:
            if gs['keep_chrM'] == 'true':
                check(kept == examined, f'{s}: fragments into bigWig (+ blacklist-removed) = pairs into MarkDuplicates',
                      f'{kept} vs {examined}')
            else:
                # keep_chrM = false drops the chrM pairs (pct_chrM of the filtered pairs, rounded to 0.01 %).
                mito = float(stats[s]['pct_chrM']) / 100 * int(stats[s]['filtered_pairs'])
                check(abs(examined - kept - mito) <= 0.0001 * examined + 1,
                      f'{s}: pairs into MarkDuplicates - fragments into bigWig = chrM pairs (keep_chrM = false)',
                      f'{examined} - {kept} = {examined - kept}, chrM ~ {mito:.0f}')
        optical = dup[s]['optical_duplicate_pairs']
        if args.test_data:
            # The test reads keep their SRA names ("SRR12246717.1"): no flow-cell coordinates.
            check(optical == 'NA', f'{s}: optical duplicates NA for SRA-named reads', optical)
        else:
            check(optical == 'NA' or optical.isdigit(), f'{s}: optical duplicate column is a count or NA', optical)

    # ---- Blacklist ----
    if blacklist_on:
        print()
        bl = {}
        with open(f'{out}/genome/blacklist.canonical.bed') as fh:
            for line in fh:
                c, b, e = line.split('\t')[:3]
                bl.setdefault(c, []).append((int(b), int(e)))
        merged_bp = sum(e - b for iv in bl.values() for b, e in iv)
        canonical = int(gs['canonical_bp'])
        check(genome_size == canonical - merged_bp, 'genome_size = canonical bp - blacklist bp',
              f'{canonical:,} - {merged_bp:,} = {genome_size:,}')
        starts = {c: [b for b, _ in sorted(iv)] for c, iv in bl.items()}
        ends = {c: [e for _, e in sorted(iv)] for c, iv in bl.items()}
        for r in sf:
            s = r['sample']
            overlap = 0
            for c, b, e in fragments(f'{out}/fragments/{s}.fragments.bed.gz'):
                if c in starts:
                    i = bisect.bisect_left(starts[c], e) - 1   # last interval starting before e
                    if i >= 0 and ends[c][i] > b:
                        overlap += 1
            check(overlap == 0, f'{s}: fragments overlapping the blacklist = 0', str(overlap))

    print(f'\n{failures} check(s) failed' if failures else '\nall checks passed')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
