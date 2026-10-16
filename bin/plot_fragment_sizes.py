#!/usr/bin/env python3
"""Plot the fragment size distributions of all samples on one figure.

Each input is <sample>.fragment_sizes.tsv with columns length, count. The y axis
is the fraction of the sample's fragments, so samples of different depth can be
compared. Writes fragment_sizes.pdf and fragment_sizes.png (0-700 bp).

Usage: plot_fragment_sizes.py a.fragment_sizes.tsv b.fragment_sizes.tsv ...
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402


def main(paths):
    fig, ax = plt.subplots(figsize=(7, 4))
    for path in sorted(paths):
        sample = os.path.basename(path).split('.fragment_sizes.tsv')[0]
        lengths, counts = [], []
        with open(path) as fh:
            next(fh)
            for line in fh:
                length, count = line.split()
                lengths.append(int(length))
                counts.append(int(count))
        total = sum(counts) or 1
        ax.plot(lengths, [c / total for c in counts], label=sample, linewidth=1)
    ax.set_xlim(0, 700)
    ax.set_xlabel('Fragment length (bp)')
    ax.set_ylabel('Fraction of fragments')
    ax.set_title('Fragment size distribution')
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig('fragment_sizes.pdf')
    fig.savefig('fragment_sizes.png', dpi=150)


if __name__ == '__main__':
    main(sys.argv[1:])
