#!/usr/bin/env python3
"""Feature sets for featureCounts and deepTools (step 11 and 14).

    make_annotations.py --gtf genes.gtf.gz --chrom-sizes canonical.chrom.sizes \
        --promoter-window 5000 [--ccre ccre.bed.gz]

Writes (only features on the canonical chromosomes):

  genes.saf          protein-coding gene bodies
  promoters.saf      protein-coding TSS +/- promoter_window, one row per gene
                     (TSS = gene start on +, gene end on -), clipped to the chromosome
  ccre.saf           ENCODE cCREs (only with --ccre)
  tss_genes.bed      protein-coding gene bodies as BED6, for computeMatrix --referencePoint TSS
  annotation_summary.tsv   row counts

SAF is featureCounts' simple annotation format: GeneID, Chr, Start, End, Strand,
1-based and end-inclusive. Protein-coding = gene_type (GENCODE) or gene_biotype
(Ensembl, dm6) "protein_coding". A gene ID seen twice (e.g. on chrX and the chrY
PAR) gets "_<chrom>" added, so every SAF row is its own feature.
"""
import argparse
import gzip
import re


def opener(path):
    with open(path, 'rb') as fh:
        gz = fh.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if gz else open(path)


def attr(field, key):
    m = re.search(rf'(?:^|;)\s*{key} "([^"]*)"', field)
    return m.group(1) if m else None


def write_saf(path, rows):
    with open(path, 'w') as fh:
        fh.write('GeneID\tChr\tStart\tEnd\tStrand\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--chrom-sizes', required=True)
    ap.add_argument('--promoter-window', type=int, required=True)
    ap.add_argument('--ccre')
    args = ap.parse_args()

    sizes = {}
    with open(args.chrom_sizes) as fh:
        for line in fh:
            c, n = line.split()[:2]
            sizes[c] = int(n)
    order = {c: i for i, c in enumerate(sorted(sizes))}

    genes, seen = [], set()
    with opener(args.gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'gene' or f[0] not in sizes:
                continue
            if (attr(f[8], 'gene_type') or attr(f[8], 'gene_biotype')) != 'protein_coding':
                continue
            gid = attr(f[8], 'gene_id')
            if gid in seen:
                gid = f'{gid}_{f[0]}'
            seen.add(gid)
            genes.append((gid, f[0], int(f[3]), min(int(f[4]), sizes[f[0]]), f[6]))
    genes.sort(key=lambda g: (order[g[1]], g[2], g[3]))
    write_saf('genes.saf', genes)

    w = args.promoter_window
    promoters = []
    for gid, c, s, e, strand in genes:
        tss = e if strand == '-' else s
        promoters.append((gid, c, max(1, tss - w), min(sizes[c], tss + w), strand))
    promoters.sort(key=lambda g: (order[g[1]], g[2], g[3]))
    write_saf('promoters.saf', promoters)

    with open('tss_genes.bed', 'w') as fh:
        for gid, c, s, e, strand in genes:
            fh.write(f'{c}\t{s - 1}\t{e}\t{gid}\t0\t{strand}\n')

    n_ccre = 'NA'
    if args.ccre:
        ccre = []
        with opener(args.ccre) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 4 and f[0] in sizes:
                    ccre.append((f[3], f[0], int(f[1]) + 1, min(int(f[2]), sizes[f[0]]), '+'))
        ccre.sort(key=lambda g: (order[g[1]], g[2], g[3]))
        write_saf('ccre.saf', ccre)
        n_ccre = len(ccre)

    with open('annotation_summary.tsv', 'w') as fh:
        fh.write(f'protein_coding_genes\t{len(genes)}\n')
        fh.write(f'promoters\t{len(promoters)}\n')
        fh.write(f'promoter_window\t{w}\n')
        fh.write(f'ccres\t{n_ccre}\n')


if __name__ == '__main__':
    main()
