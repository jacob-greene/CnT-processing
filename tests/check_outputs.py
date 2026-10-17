#!/usr/bin/env python3
"""Known-answer checks on a finished CnT-processing run (steps 1-16).

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
    total pairs = MANIFEST.tsv (test data; IgG_rep1 has two lanes),
    filtered <= aligned <= after trimming <= total
  * duplicates kept: fragments (+ blacklist-removed) = pairs examined by MarkDuplicates
    (--test-data only: every test pair is on chr22/chrM with a fragment < 1000 bp)
  * optical duplicates are NA for the SRA-named test reads (--test-data)
  * blacklist_filter = true: 0 fragments overlap the blacklist, and genome_size =
    canonical bp - merged blacklist bp
  * featureCounts: every matrix value = the number of fragments overlapping that
    feature, counted independently from the fragment BED and the SAF; Assigned =
    fragments overlapping >= 1 feature; summary total = fragments; rate > 0;
    promoter rows = protein-coding genes in the GTF; promoter = TSS +/- window
  * peaks_frip.tsv: peak count and widths from the SEACR BED; FRiP and
    control_frip recomputed from the fragment BEDs; SEACR mode = control for
    every target with a control; control_frip < FRiP (test data)
  * deepTools: the peak and TSS plots exist; the TSS matrix has one row per
    protein-coding gene
  * correlation: the 5 kb bin means rebuilt from the fragment BEDs x s_cov equal
    the pipeline's bin file and the bigWigs; the matrix = Pearson of log1p of these
    bins (bins that are 0 in every sample skipped)
  * MultiQC: the report exists and holds every expected section
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


# ---------------------------------------------------------------------------
# Steps 11-16
# ---------------------------------------------------------------------------

def frag_arrays(path):
    """Fragments per chromosome as sorted numpy arrays of starts and of ends (0-based, half-open)."""
    by = {}
    for c, b, e in fragments(path):
        by.setdefault(c, ([], []))
        by[c][0].append(b)
        by[c][1].append(e)
    return {c: (np.sort(np.array(b)), np.sort(np.array(e))) for c, (b, e) in by.items()}


def per_feature_counts(fr, feats):
    """feats: list of (chrom, start0, end). Fragments [s, e) overlapping [S, E) by >= 1 bp:
    s < E and e > S = (#starts < E) - (#ends <= S), as no fragment ends before it starts."""
    out = np.zeros(len(feats), dtype=np.int64)
    for i, (c, S, E) in enumerate(feats):
        if c in fr:
            st, en = fr[c]
            out[i] = np.searchsorted(st, E, 'left') - np.searchsorted(en, S, 'right')
    return out


def bin_sums(starts, ends, values, width, n_bins):
    """Sum of value x overlap bp per bin [k*width, (k+1)*width) for intervals [start, end)."""
    k0, k1 = starts // width, (ends - 1) // width
    out = np.zeros(n_bins + 1)
    one = k0 == k1
    np.add.at(out, k0[one], values[one] * (ends[one] - starts[one]))
    m = ~one
    np.add.at(out, k0[m], values[m] * ((k0[m] + 1) * width - starts[m]))
    np.add.at(out, k1[m], values[m] * (ends[m] - k1[m] * width))
    full = np.zeros(n_bins + 2)                      # bins strictly between k0 and k1
    np.add.at(full, k0[m] + 1, values[m])
    np.add.at(full, k1[m], -values[m])
    return (out + np.cumsum(full)[:n_bins + 1] * width)[:n_bins]


def merge(intervals):
    """Sorted, merged (start, end) arrays per chromosome."""
    by = {}
    for c, b, e in intervals:
        by.setdefault(c, []).append((b, e))
    res = {}
    for c, iv in by.items():
        iv.sort()
        ms, me = [], []
        for b, e in iv:
            if ms and b <= me[-1]:
                me[-1] = max(me[-1], e)
            else:
                ms.append(b)
                me.append(e)
        res[c] = (np.array(ms), np.array(me))
    return res


def n_overlapping(path, merged):
    """Fragments in a BED overlapping >= 1 merged interval by >= 1 bp, and all fragments."""
    hit = total = 0
    by = {}
    for c, b, e in fragments(path):
        by.setdefault(c, ([], []))
        by[c][0].append(b)
        by[c][1].append(e)
        total += 1
    for c, (b, e) in by.items():
        if c not in merged:
            continue
        ms, me = merged[c]
        b, e = np.array(b), np.array(e)
        i = np.searchsorted(me, b, 'right')          # first merged interval ending after the fragment start
        ok = i < len(ms)
        hit += int((ms[np.minimum(i, len(ms) - 1)] < e)[ok].sum())
    return hit, total


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3:
                rows.append((f[0], int(f[1]), int(f[2])))
    return rows


def protein_coding_genes(gtf, sizes):
    n = 0
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) > 8 and f[2] == 'gene' and f[0] in sizes and \
                    ('gene_type "protein_coding"' in f[8] or 'gene_biotype "protein_coding"' in f[8]):
                n += 1
    return n


def downstream_checks(out, sf, sizes, test_data):
    samples = [r['sample'] for r in sf]
    frag_path = {s: f'{out}/fragments/{s}.fragments.bed.gz' for s in samples}
    n_frag = {r['sample']: int(r['fragments']) for r in sf}

    # ---- featureCounts ----
    print()
    fc_dir = f'{out}/featurecounts'
    sets = [x for x in ('genes', 'promoters', 'ccre') if os.path.exists(f'{fc_dir}/{x}.counts.tsv')]
    check('genes' in sets and 'promoters' in sets, 'featureCounts: gene and promoter matrices exist', ', '.join(sets))
    fr = {s: frag_arrays(frag_path[s]) for s in samples}
    for name in sets:
        saf = read_tsv(f'{fc_dir}/annotation/{name}.saf')
        feats = [(r['Chr'], int(r['Start']) - 1, int(r['End'])) for r in saf]
        mat = read_tsv(f'{fc_dir}/{name}.counts.tsv')
        check([r['Geneid'] for r in mat] == [r['GeneID'] for r in saf],
              f'featureCounts {name}: matrix rows = SAF rows', f'{len(mat):,} features')
        summary = {}
        with open(f'{fc_dir}/{name}.counts.tsv.summary') as fh:
            head = next(fh).rstrip('\n').split('\t')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                summary[f[0]] = dict(zip(head[1:], map(int, f[1:])))
        merged = merge(feats)
        for s in samples:
            want = per_feature_counts(fr[s], feats)
            got = np.array([int(r[s]) for r in mat])
            check(np.array_equal(got, want), f'featureCounts {name} {s}: every count = fragments overlapping the feature',
                  f'{int((got != want).sum())} of {len(got):,} differ; column sum {got.sum():,}')
            total = sum(v[s] for v in summary.values())
            hit, _ = n_overlapping(frag_path[s], merged)
            assigned = summary['Assigned'][s]
            check(total == n_frag[s], f'featureCounts {name} {s}: summary total = fragments', f'{total:,}')
            check(assigned == hit and assigned > 0,
                  f'featureCounts {name} {s}: Assigned = fragments overlapping >= 1 feature (> 0)',
                  f'{assigned:,} vs {hit:,}, rate {assigned / max(1, total):.3f}')

    ann = read_kv(f'{fc_dir}/annotation/annotation_summary.tsv')
    n_pc = protein_coding_genes(f'{out}/genome/genes.gtf.gz', sizes)
    prom = read_tsv(f'{fc_dir}/annotation/promoters.saf')
    genes = {r['GeneID']: r for r in read_tsv(f'{fc_dir}/annotation/genes.saf')}
    check(len(prom) == n_pc, 'promoter SAF rows = protein-coding genes (TSSs) in the GTF', f'{len(prom):,} vs {n_pc:,}')
    w = int(ann['promoter_window'])
    bad = 0
    for p in prom:
        g = genes[p['GeneID']]
        tss = int(g['End']) if g['Strand'] == '-' else int(g['Start'])
        if (int(p['Start']), int(p['End'])) != (max(1, tss - w), min(sizes[p['Chr']], tss + w)):
            bad += 1
    check(bad == 0, f'every promoter = its gene TSS +/- {w} bp (clipped to the chromosome)', f'{bad} wrong')

    # ---- SEACR peaks + peaks_frip.tsv ----
    print()
    pf = read_tsv(f'{out}/peaks_frip.tsv')
    controls = {r['control'] for r in pf if r['control'] != 'NA'}
    check(not (controls & {r['sample'] for r in pf}), 'controls are not peak-called (no peaks_frip row)')
    if test_data:
        check(len(pf) == 2 and all(r['seacr_mode'] == 'control' for r in pf),
              'test data: both targets peak-called in control mode', ', '.join(f'{r["sample"]}={r["seacr_mode"]}' for r in pf))
    for r in pf:
        s = r['sample']
        has_ctrl = r['control'] != 'NA'
        check(r['seacr_mode'] == ('control' if has_ctrl else 'top0.01'), f'{s}: SEACR mode matches the control column',
              r['seacr_mode'])
        peak_sets = [('', f'{out}/seacr/{s}.stringent.bed')]
        if has_ctrl:
            peak_sets.append(('top01_', f'{out}/seacr/top0.01/{s}.stringent.bed'))
        for prefix, path in peak_sets:
            peaks = read_bed(path)
            widths = sorted(e - b for _, b, e in peaks)
            label = f'{s} {"top0.01" if prefix else "primary"} peaks'
            check(int(r[prefix + 'peaks']) == len(peaks), f'{label}: peak count = SEACR BED lines', f'{len(peaks):,}')
            if widths:
                check(float(r[prefix + 'median_width']) == float(np.median(widths)) and
                      abs(float(r[prefix + 'mean_width']) - np.mean(widths)) < 0.051 and
                      int(r[prefix + 'bp_in_peaks']) == sum(widths),
                      f'{label}: median / mean width and bp in peaks', f'{np.median(widths):g} / {np.mean(widths):.1f} / {sum(widths):,}')
            merged = merge(peaks)
            hit, total = n_overlapping(frag_path[s], merged)
            frip = hit / total
            check(abs(float(r[prefix + 'frip']) - frip) < 1e-6, f'{label}: FRiP recomputed from the fragment BED', f'{frip:.4f}')
            if has_ctrl:
                chit, ctot = n_overlapping(frag_path[r['control']], merged)
                cfrip = chit / ctot
                check(abs(float(r[prefix + 'control_frip']) - cfrip) < 1e-6,
                      f'{label}: control_frip recomputed ({r["control"]} fragments)', f'{cfrip:.4f}')
                if not prefix and test_data:
                    # Not a known answer here: the test IgG is ~88% chrM pairs (pairs were chosen
                    # on chr22 + chrM), so a target peak on chrM holds most IgG fragments.
                    print(f'INFO  {s}: control_frip {cfrip:.4f}, FRiP {frip:.4f} (test data: not checked)')
                elif not prefix:
                    check(cfrip < frip, f'{s}: control_frip < FRiP', f'{cfrip:.4f} < {frip:.4f}')

    # ---- deepTools ----
    print()
    dt = f'{out}/deeptools'
    for r in pf:
        s = r['sample']
        if int(r['peaks']) > 0:
            ok = all(os.path.getsize(f'{dt}/peaks/{s}.peaks_{k}.pdf') > 0 for k in ('heatmap', 'profile'))
            check(ok, f'{s}: peak heatmap and profile written')
    kinds = ['coverage_norm'] + (['readcount_scaled'] if os.path.isdir(f'{out}/bigwig/readcount_scaled') else [])
    for k in kinds:
        check(os.path.getsize(f'{dt}/tss/tss_profile_{k}.pdf') > 0, f'TSS profile written: tss_profile_{k}.pdf')
    check(os.path.getsize(f'{dt}/tss/tss_heatmap_coverage_norm.pdf') > 0, 'TSS heatmap written')
    with gzip.open(f'{dt}/tss/tss_matrix_coverage_norm.gz', 'rt') as fh:
        import json
        head = json.loads(next(fh)[1:])
        rows = sum(1 for _ in fh)
    n_tss = sum(1 for _ in open(f'{fc_dir}/annotation/tss_genes.bed'))
    check(rows == n_tss and head['sample_labels'] == samples,
          'TSS matrix: one row per protein-coding gene, every sample a column', f'{rows:,} rows vs {n_tss:,} genes')

    # ---- Correlation ----
    # The bin values are recomputed from the fragment BEDs (true 5 kb means x s_cov) and
    # compared with the pipeline's bin file and with the bigWigs, so a wrong bigWig or a
    # wrong bin file both turn this red. The matrix is then recomputed from these bins.
    print()
    with gzip.open(f'{out}/correlation/correlation_bins.tsv.gz', 'rt') as fh:
        labels = [x.strip("'") for x in next(fh).rstrip('\n').split('\t')[3:]]
        rows = [line.rstrip('\n').split('\t') for line in fh]
    chroms = np.array([r[0] for r in rows])
    b_start = np.array([int(r[1]) for r in rows])
    b_end = np.array([int(r[2]) for r in rows])
    pipe = np.nan_to_num(np.array([[float(x) for x in r[3:]] for r in rows]))   # NaN = no data = 0
    width = int(np.max(b_end - b_start))
    factor = {r['sample']: float(r['s_cov']) for r in sf}
    want_bins = np.zeros_like(pipe)
    bw_bins = np.zeros_like(pipe)
    for j, s in enumerate(labels):
        fr = frag_arrays(f'{out}/fragments/{s}.fragments.bed.gz')
        bw = pyBigWig.open(f'{out}/bigwig/coverage_norm/{s}.coverage_norm.bw')
        for c in np.unique(chroms):
            idx = np.flatnonzero(chroms == c)
            k, n_bins = b_start[idx] // width, int(b_start[idx].max() // width) + 1
            if c in fr:
                st, en = fr[c]
                want_bins[idx, j] = bin_sums(st, en, np.ones(len(st)), width, n_bins)[k]
            iv = bw.intervals(c) or []
            if iv:
                a = np.array(iv, dtype=float)
                bw_bins[idx, j] = bin_sums(a[:, 0].astype(np.int64), a[:, 1].astype(np.int64), a[:, 2],
                                           width, n_bins)[k]
        bw.close()
        length = (b_end - b_start).astype(float)
        want_bins[:, j] *= factor[s] / length
        bw_bins[:, j] /= length
    for j, s in enumerate(labels):
        for name, got_bins in (('pipeline bin file', pipe), ('bigWig', bw_bins)):
            err = np.abs(got_bins[:, j] - want_bins[:, j]) / np.maximum(want_bins[:, j], 1)
            check(err.max() < 1e-4, f'{s}: 5 kb bin means from the {name} = fragment BED x s_cov',
                  f'{len(err):,} bins, max error {err.max():.1e}')
    vals = want_bins[(want_bins != 0).any(axis=1)]   # --skipZeros: drop bins that are 0 in every sample
    want = np.corrcoef(np.log1p(vals).T)
    got = {}
    with open(f'{out}/correlation/correlation_matrix.tsv') as fh:
        lines = [l.rstrip('\n').split('\t') for l in fh if not l.startswith('#')]
    cols = [x.strip("'") for x in lines[0][1:]]
    for f in lines[1:]:
        for c, v in zip(cols, f[1:]):
            got[(f[0].strip("'"), c)] = float(v)
    err = max(abs(got[(a, b)] - want[i, j]) for i, a in enumerate(labels) for j, b in enumerate(labels))
    check(sorted(labels) == sorted(samples) and err < 1e-4,
          'correlation matrix = Pearson of log1p(5 kb bin means), bins rebuilt from the fragment BEDs',
          f'{len(vals):,} bins, max error {err:.1e}')
    for i, a in enumerate(labels):
        best = max((want[i, j], b) for j, b in enumerate(labels) if j != i)
        print(f'INFO  {a}: most correlated with {best[1]} (r = {best[0]:.3f})')

    # ---- MultiQC ----
    print()
    mq = f'{out}/multiqc'
    check(os.path.getsize(f'{mq}/multiqc_report.html') > 0, 'MultiQC report written')
    html = open(f'{mq}/multiqc_report.html', encoding='utf-8', errors='replace').read().lower()  # anchors vary in case
    want_ids = ['fastqc', 'cutadapt', 'bowtie2', 'bowtie2_ecoli', 'samtools', 'picard', 'featurecounts_genes',
                'featurecounts_promoters', 'alignment_stats', 'duplication', 'scale_factors', 'peaks_frip',
                'fragment_sizes', 'correlation_heatmap'] + (['featurecounts_ccre'] if 'ccre' in sets else [])
    missing = [x for x in want_ids if f'id="{x}"' not in html and f"id='{x}'" not in html]
    check(not missing, 'MultiQC report holds every expected section', 'missing: ' + ', '.join(missing) if missing else
          f'{len(want_ids)} sections')


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
    manifest = read_kv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'MANIFEST.tsv'))
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
            # IgG_rep1 has two lanes with the same file names in different folders:
            # this holds only if both lanes were joined.
            want = manifest.get(f'{s}.ci_pairs')
            check(want is not None and tot == int(want), f'{s}: total pairs = tests/data/MANIFEST.tsv',
                  f'{tot} vs {want} ({manifest.get(f"{s}.ci_lanes", "1 lane")})')

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

    downstream_checks(out, sf, sizes, args.test_data)

    print(f'\n{failures} check(s) failed' if failures else '\nall checks passed')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
