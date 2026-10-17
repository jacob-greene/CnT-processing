#!/usr/bin/env python3
"""Negative control for the per-base check in check_outputs.py.

    python3 tests/mutate_bigwig.py <in.bw> <out.bw>

Shuffles the values among bigWig intervals of the same length on the same
chromosome. Every value moves, but the total signal stays exactly the same,
so a check of totals still passes and the per-base check must fail.
check_outputs.py imports shuffle() and runs this control on every call.
"""
import sys

import numpy as np
import pyBigWig


def shuffle(src, dst, seed=7):
    """Write dst = src with values permuted among equal-length intervals. Returns # values moved."""
    bw = pyBigWig.open(src)
    chroms = bw.chroms()
    out = pyBigWig.open(dst, 'w')
    out.addHeader(list(chroms.items()))
    rng = np.random.default_rng(seed)
    moved = 0
    for c in chroms:
        iv = bw.intervals(c) or []
        vals = np.array([v for _, _, v in iv])
        lengths = np.array([e - b for b, e, _ in iv])
        for n in np.unique(lengths):
            idx = np.flatnonzero(lengths == n)
            new = vals[idx][rng.permutation(len(idx))]
            moved += int((new != vals[idx]).sum())
            vals[idx] = new
        if iv:
            out.addEntries([c] * len(iv), [b for b, _, _ in iv], ends=[e for _, e, _ in iv],
                           values=[float(v) for v in vals])
    bw.close()
    out.close()
    return moved


if __name__ == '__main__':
    src, dst = sys.argv[1:3]
    print(f'{src}: {shuffle(src, dst)} interval values moved, total unchanged')
