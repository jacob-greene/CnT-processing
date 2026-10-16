#!/usr/bin/env bash
# Compute the canonical chromosome sizes and the effective genome size.
#
# Usage: genome_size.sh <chrom.sizes> <mito_contig> <keep_chrM true|false> [blacklist.bed[.gz]]
#
# Canonical chromosomes = every name in the UCSC chrom.sizes file that has no "_"
# (this drops _random, chrUn_*, _alt, _fix and _hap contigs). The mitochondrial
# contig is kept or dropped by <keep_chrM>.
#
#   genome_size = sum(canonical chromosome sizes)                    (no blacklist)
#   genome_size = sum(canonical sizes) - merged blacklist bp         (with blacklist)
#
# Blacklist intervals are clipped to the canonical chromosomes and merged first, so
# overlapping intervals are not counted twice.
#
# Writes: canonical.chrom.sizes (sorted by name), genome_size.tsv (key<TAB>value),
#         and blacklist.canonical.bed when a blacklist is given.
set -euo pipefail

sizes=$1
mito=$2
keep_mito=$3
blacklist=${4:-}

awk -v OFS='\t' -v mito="$mito" -v keep="$keep_mito" '
    $1 ~ /_/ { next }
    $1 == mito && keep != "true" { next }
    { print $1, $2 }
' "$sizes" | sort -k1,1 > canonical.chrom.sizes

canonical_bp=$(awk '{ s += $2 } END { printf "%d", s }' canonical.chrom.sizes)
mito_bp=$(awk -v mito="$mito" '$1 == mito { print $2 }' canonical.chrom.sizes)

{
    printf 'chromosomes\t%s\n' "$(wc -l < canonical.chrom.sizes | tr -d ' ')"
    printf 'mito_contig\t%s\n' "$mito"
    printf 'keep_chrM\t%s\n' "$keep_mito"
    printf 'mito_bp_included\t%s\n' "${mito_bp:-0}"
    printf 'canonical_bp\t%s\n' "$canonical_bp"
} > genome_size.tsv

if [[ -n "$blacklist" ]]; then
    # Clip to canonical chromosomes, then merge overlapping intervals.
    zcat -f "$blacklist" \
        | awk -v OFS='\t' 'NR == FNR { len[$1] = $2; next }
                           ($1 in len) { s = ($2 < 0 ? 0 : $2); e = ($3 > len[$1] ? len[$1] : $3);
                                         if (e > s) print $1, s, e }' canonical.chrom.sizes - \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - > blacklist.canonical.bed
    blacklist_bp=$(awk '{ s += $3 - $2 } END { printf "%d", s }' blacklist.canonical.bed)
    printf 'blacklist_filter\ttrue\n'                 >> genome_size.tsv
    printf 'blacklist_bp\t%s\n' "$blacklist_bp"        >> genome_size.tsv
    printf 'genome_size\t%s\n' "$((canonical_bp - blacklist_bp))" >> genome_size.tsv
else
    printf 'blacklist_filter\tfalse\n'                >> genome_size.tsv
    printf 'genome_size\t%s\n' "$canonical_bp"         >> genome_size.tsv
fi

cat genome_size.tsv
