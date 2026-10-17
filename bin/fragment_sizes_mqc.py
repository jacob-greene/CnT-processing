#!/usr/bin/env python3
"""Fragment size tables -> one MultiQC custom-content line graph (JSON on stdout).

    fragment_sizes_mqc.py <sample>[.fragment_sizes].tsv ...   (columns: length, count)
"""
import json
import os
import sys

data = {}
for path in sys.argv[1:]:
    sample = os.path.basename(path)
    for ext in ('.fragment_sizes.tsv', '.tsv'):
        if sample.endswith(ext):
            sample = sample[: -len(ext)]
            break
    with open(path) as fh:
        next(fh)
        data[sample] = {int(n): int(c) for n, c in (line.split('\t')[:2] for line in fh) if int(n) <= 700}
json.dump({
    'id': 'fragment_sizes',
    'section_name': 'Fragment sizes',
    'description': 'Fragment length (|TLEN| of properly paired reads after filtering), 0-700 bp.',
    'plot_type': 'linegraph',
    'pconfig': {'id': 'fragment_sizes_plot', 'title': 'Fragment sizes',
                'xlab': 'Fragment length (bp)', 'ylab': 'Read pairs', 'xmin': 0, 'xmax': 700},
    'data': data,
}, sys.stdout)
