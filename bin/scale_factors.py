#!/usr/bin/env python3
"""Write scale_factors.tsv: the bigWig scale factors of every sample.

    C_i     = total per-base coverage of sample i = sum of its fragment lengths
    s_cov_i = genome_size / C_i                      (coverage-normalized bigWig)
    s_rc_i  = genome_size / mean_j(C_j)              (read-count scaled bigWig;
                                                      j = every sample in i's pool)
    ratio_i = C_i / mean_j(C_j)                      (= s_rc_i / s_cov_i)

Usage:
    scale_factors.py --genome-size genome_size.tsv --blacklist-filter true|false \
        --pools id1=pool1,id2=pool1 id1.fragments.tsv id2.fragments.tsv ...
"""
import argparse
import os
import sys


def read_kv(path):
    with open(path) as fh:
        return dict(line.rstrip('\n').split('\t', 1) for line in fh if '\t' in line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--genome-size', required=True)
    ap.add_argument('--blacklist-filter', required=True, choices=['true', 'false'])
    ap.add_argument('--pools', required=True)
    ap.add_argument('stats', nargs='+')
    args = ap.parse_args()

    genome_size = int(read_kv(args.genome_size)['genome_size'])
    pools = dict(p.split('=', 1) for p in args.pools.split(','))

    samples = {}
    for path in args.stats:
        sample = os.path.basename(path)[: -len('.fragments.tsv')]
        samples[sample] = read_kv(path)
    missing = set(pools) - set(samples)
    if missing:
        sys.exit(f'ERROR: no fragment stats for: {", ".join(sorted(missing))}')

    coverage = {s: int(samples[s]['coverage_bp']) for s in pools}
    for s, c in coverage.items():
        if c == 0:
            sys.exit(f'ERROR: sample {s} has no fragments left; cannot scale its bigWig')

    pool_mean = {}
    for pool in set(pools.values()):
        members = [s for s in pools if pools[s] == pool]
        pool_mean[pool] = sum(coverage[s] for s in members) / len(members)

    cols = ['sample', 'pool', 'fragments', 'C_i', 's_cov', 's_rc', 'C_i_over_pool_mean',
            'blacklist_filter', 'fragments_removed_blacklist', 'genome_size']
    print('\t'.join(cols))
    for s in pools:  # sample sheet order
        mean = pool_mean[pools[s]]
        row = [s, pools[s], samples[s]['fragments'], str(coverage[s]),
               f'{genome_size / coverage[s]:.10g}', f'{genome_size / mean:.10g}',
               f'{coverage[s] / mean:.6f}', args.blacklist_filter,
               samples[s]['fragments_removed_blacklist'], str(genome_size)]
        print('\t'.join(row))


if __name__ == '__main__':
    main()
