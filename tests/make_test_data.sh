#!/usr/bin/env bash
# Build the test data for CnT-processing.
#
# Source: K562 CUT&Tag from Kaya-Okur et al. 2020 (GEO GSE145187), the samples
# used by the Henikoff lab CUT&Tag tutorial. One SRA run per sample:
#
#   sample         SRA experiment  SRA run      role
#   H3K27me3_rep1  SRX8754646      SRR12246717  target
#   H3K4me3_rep1   SRX7713692      SRR11074254  target
#   IgG_rep1       SRX8468909      SRR11923224  control (is_control=true)
#
# Two tiers are built:
#
#   smoke  ~1M read pairs per sample (seqtk sample -s 100), for the full-hg38
#          run on the cluster. Written to $WORK/smoke/ (not committed).
#
#   CI     the smoke reads restricted to pairs that align concordantly to chr22
#          or chrM in hg38, plus pairs that align concordantly to E. coli
#          (so the carry-over QC is non-zero). Each sample is then downsampled
#          to a different depth (ratio 4 : 2 : 1) so read-count scaling is tested.
#          Written to $OUT/ (committed, except the FASTA and Bowtie2 index):
#            fastq/<sample>_R{1,2}.fastq.gz
#            samplesheet.csv
#            genome/chr22_chrM.chrom.sizes
#            genome/chr22_chrM.gencode.v50.gtf.gz
#            genome/chr22_chrM.ccre.bed.gz
#            genome/chr22_chrM.blacklist.v2.bed.gz
#            genome/chr22_chrM.fa.gz          (UCSC chr22 + chrM, gzip -n -9)
#            bowtie2/chr22_chrM.*.bt2         (git-ignored)
#            MANIFEST.tsv                     (sources, checksums, read counts)
#
# Requirements on PATH: curl, md5sum, gzip, awk, bowtie2, bowtie2-build,
# samtools, seqtk. At Fred Hutch, load them in this order (later modules
# replace the GCC toolchain that seqtk loads):
#
#   ml seqtk/1.3-GCC-10.2.0
#   ml Bowtie2/2.5.4-GCC-13.2.0 SAMtools/1.19.2-GCC-13.2.0
#
# Usage:
#   HG38_BT2=/path/to/hg38/bowtie2/prefix bash tests/make_test_data.sh
#
# Environment variables:
#   HG38_BT2   full hg38 (UCSC names) Bowtie2 index prefix, used only to select
#              chr22/chrM pairs. Default: the Fred Hutch shared iGenomes index.
#   WORK       scratch folder for downloads and intermediate files
#              (default: ./test_data_work)
#   OUT        output folder for the CI data (default: <repo>/tests/data)
#   THREADS    threads for Bowtie2 (default: 8)
#   SMOKE_PAIRS  read pairs per sample in the smoke tier (default: 1000000)
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
HG38_BT2=${HG38_BT2:-/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome}
WORK=${WORK:-$PWD/test_data_work}
OUT=${OUT:-$REPO/tests/data}
THREADS=${THREADS:-8}
SMOKE_PAIRS=${SMOKE_PAIRS:-1000000}
SEED=100

# Pinned sources.
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38
GTF_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
CCRE_URL=https://www.encodeproject.org/files/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz
BLACKLIST_URL=https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz
ECOLI_URL="https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true"

# sample  run  ENA path  R1 md5  R2 md5   (from the ENA filereport API)
SAMPLES=(
  "H3K27me3_rep1 SRR12246717 SRR122/017/SRR12246717 cebaa8e610ce4e325956c02f555ff41b 14110fe972c02c592ff8f24eee1268aa"
  "H3K4me3_rep1  SRR11074254 SRR110/054/SRR11074254 9eef41a94fa6f63973f1160e8aeee9fe 429e80630533759e6b3f09075fbcd5d9"
  "IgG_rep1      SRR11923224 SRR119/024/SRR11923224 af973d0fb466680bde1deda33befeeb0 be7ea67f8e5abd9bba3541796d90f27f"
)
# CI depth of each sample, as a fraction of the smallest chr22+chrM+E. coli pair count.
CI_FRACTION=(1 0.5 0.25)
CONTIGS="chr22 chrM"

BT2_OPTS="--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700"

for tool in curl md5sum gzip awk bowtie2 bowtie2-build samtools seqtk; do
    command -v "$tool" > /dev/null || { echo "ERROR: '$tool' is not on PATH" >&2; exit 1; }
done
[[ -e "$HG38_BT2.1.bt2" || -e "$HG38_BT2.1.bt2l" ]] || {
    echo "ERROR: no Bowtie2 index at HG38_BT2=$HG38_BT2" >&2; exit 1; }

mkdir -p "$WORK"/{raw,smoke,align,ref} "$OUT"/{fastq,genome,bowtie2}
log() { echo "[$(date +%H:%M:%S)] $*"; }

# fetch <url> <dest> [md5]
fetch() {
    if [[ ! -s "$2" ]]; then
        log "download $1"
        curl -sfL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
    fi
    if [[ -n "${3:-}" ]]; then
        echo "$3  $2" | md5sum -c --quiet -
    fi
}

# count_pairs <R1.fastq.gz>
count_pairs() { echo $(( $(zcat "$1" | wc -l) / 4 )); }

# ---------------------------------------------------------------------------
# 1. Raw FASTQs from ENA, checked against the ENA md5 sums
# ---------------------------------------------------------------------------
for entry in "${SAMPLES[@]}"; do
    read -r sample run path md5_1 md5_2 <<< "$entry"
    fetch "$ENA/${path}_1.fastq.gz" "$WORK/raw/${run}_1.fastq.gz" "$md5_1"
    fetch "$ENA/${path}_2.fastq.gz" "$WORK/raw/${run}_2.fastq.gz" "$md5_2"
done

# ---------------------------------------------------------------------------
# 2. Smoke tier: ~1M pairs per sample. Same seed and count on R1 and R2 keeps
#    the mates paired.
# ---------------------------------------------------------------------------
for entry in "${SAMPLES[@]}"; do
    read -r sample run _ <<< "$entry"
    for r in 1 2; do
        dest="$WORK/smoke/${sample}_R$r.fastq.gz"
        [[ -s "$dest" ]] || seqtk sample -s "$SEED" "$WORK/raw/${run}_$r.fastq.gz" "$SMOKE_PAIRS" | gzip -c > "$dest"
    done
done
{
    echo "sample,fastq_1,fastq_2,control,is_control,pool"
    echo "H3K27me3_rep1,H3K27me3_rep1_R1.fastq.gz,H3K27me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1"
    echo "H3K4me3_rep1,H3K4me3_rep1_R1.fastq.gz,H3K4me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1"
    echo "IgG_rep1,IgG_rep1_R1.fastq.gz,IgG_rep1_R2.fastq.gz,,true,pool1"
} > "$WORK/smoke/samplesheet.csv"

# ---------------------------------------------------------------------------
# 3. Test references: chr22 + chrM FASTA, chrom sizes, GTF, cCREs, blacklist
# ---------------------------------------------------------------------------
fetch "$UCSC/bigZips/hg38.chrom.sizes" "$WORK/ref/hg38.chrom.sizes"
fetch "$GTF_URL"       "$WORK/ref/gencode.v50.annotation.gtf.gz"
fetch "$CCRE_URL"      "$WORK/ref/ENCFF420VPZ.bed.gz"      ea1820bd71581ba01c94fc600fbf794d
fetch "$BLACKLIST_URL" "$WORK/ref/hg38-blacklist.v2.bed.gz" 83fe6bf8187a64dee8079b80f75ba289
fetch "$ECOLI_URL"     "$WORK/ref/U00096.3.fa"
for c in $CONTIGS; do
    fetch "$UCSC/chromosomes/$c.fa.gz" "$WORK/ref/$c.fa.gz"
done

G="$OUT/genome/chr22_chrM"
fetch "$UCSC/chromosomes/md5sum.txt" "$WORK/ref/chromosomes.md5sum.txt"
for c in $CONTIGS; do
    grep " $c.fa.gz\$" "$WORK/ref/chromosomes.md5sum.txt" | sed "s|$c.fa.gz|$WORK/ref/$c.fa.gz|" | md5sum -c --quiet -
done
# -n: no file name or time stamp, so a rebuild gives the same bytes.
zcat $(for c in $CONTIGS; do echo "$WORK/ref/$c.fa.gz"; done) | gzip -n -9 -c > "$G.fa.gz"
awk -v keep="$CONTIGS" 'BEGIN { n = split(keep, k, " "); for (i = 1; i <= n; i++) want[k[i]] = 1 }
                        ($1 in want)' "$WORK/ref/hg38.chrom.sizes" | sort -k1,1 > "$G.chrom.sizes"
# Keep the GTF header lines and every feature on the test contigs.
zcat "$WORK/ref/gencode.v50.annotation.gtf.gz" \
    | awk -F '\t' '/^#/ || $1 == "chr22" || $1 == "chrM"' | gzip -c > "$G.gencode.v50.gtf.gz"
zcat "$WORK/ref/ENCFF420VPZ.bed.gz" \
    | awk -F '\t' '$1 == "chr22" || $1 == "chrM"' | gzip -c > "$G.ccre.bed.gz"
zcat "$WORK/ref/hg38-blacklist.v2.bed.gz" \
    | awk -F '\t' '$1 == "chr22" || $1 == "chrM"' | gzip -c > "$G.blacklist.v2.bed.gz"

[[ -s "$OUT/bowtie2/chr22_chrM.1.bt2" ]] || \
    bowtie2-build --threads "$THREADS" -q "$G.fa.gz" "$OUT/bowtie2/chr22_chrM"
[[ -s "$WORK/ref/ecoli.1.bt2" ]] || \
    bowtie2-build --threads "$THREADS" -q "$WORK/ref/U00096.3.fa" "$WORK/ref/ecoli"

# ---------------------------------------------------------------------------
# 4. CI tier: select pairs that are concordant on chr22/chrM (full hg38 index)
#    or concordant on E. coli, then downsample to unequal depths.
# ---------------------------------------------------------------------------
declare -A available
for entry in "${SAMPLES[@]}"; do
    read -r sample _ <<< "$entry"
    names="$WORK/align/$sample.keep.txt"
    if [[ ! -s "$names" ]]; then
        log "align $sample to hg38 and E. coli"
        R1="$WORK/smoke/${sample}_R1.fastq.gz"; R2="$WORK/smoke/${sample}_R2.fastq.gz"
        # -f 2: properly paired; -F 0x904: drop unmapped, secondary, supplementary.
        bowtie2 $BT2_OPTS -p "$THREADS" -x "$HG38_BT2" -1 "$R1" -2 "$R2" 2> "$WORK/align/$sample.hg38.log" \
            | samtools view -f 2 -F 0x904 - \
            | awk -F '\t' '($3 == "chr22" || $3 == "chrM") && $7 == "=" { print $1 }' \
            | sort -u > "$WORK/align/$sample.hg38.names"
        bowtie2 $BT2_OPTS --no-overlap --no-dovetail -p "$THREADS" -x "$WORK/ref/ecoli" -1 "$R1" -2 "$R2" \
                2> "$WORK/align/$sample.ecoli.log" \
            | samtools view -f 2 -F 0x904 - | cut -f 1 | sort -u > "$WORK/align/$sample.ecoli.names"
        cat "$WORK/align/$sample.hg38.names" "$WORK/align/$sample.ecoli.names" | sort -u > "$names"
    fi
    available[$sample]=$(wc -l < "$names")
done

min_pairs=$(for s in "${!available[@]}"; do echo "${available[$s]}"; done | sort -n | head -1)
log "pairs available (chr22 + chrM + E. coli): $(for s in "${!available[@]}"; do printf '%s=%s ' "$s" "${available[$s]}"; done)"

i=0
for entry in "${SAMPLES[@]}"; do
    read -r sample _ <<< "$entry"
    target=$(awk -v n="$min_pairs" -v f="${CI_FRACTION[$i]}" 'BEGIN { printf "%d", n * f }')
    for r in 1 2; do
        # seqtk subseq keeps file order, so R1 and R2 stay in step; the same
        # seed and count in seqtk sample then keeps the mates paired.
        seqtk subseq "$WORK/smoke/${sample}_R$r.fastq.gz" "$WORK/align/$sample.keep.txt" \
            | seqtk sample -s "$SEED" - "$target" | gzip -n -c > "$OUT/fastq/${sample}_R$r.fastq.gz"
    done
    i=$((i + 1))
done

{
    echo "sample,fastq_1,fastq_2,control,is_control,pool"
    echo "H3K27me3_rep1,fastq/H3K27me3_rep1_R1.fastq.gz,fastq/H3K27me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1"
    echo "H3K4me3_rep1,fastq/H3K4me3_rep1_R1.fastq.gz,fastq/H3K4me3_rep1_R2.fastq.gz,IgG_rep1,false,pool1"
    echo "IgG_rep1,fastq/IgG_rep1_R1.fastq.gz,fastq/IgG_rep1_R2.fastq.gz,,true,pool1"
} > "$OUT/samplesheet.csv"

# ---------------------------------------------------------------------------
# 5. Manifest: what was built, from what, and how many reads
# ---------------------------------------------------------------------------
{
    printf 'item\tvalue\n'
    for entry in "${SAMPLES[@]}"; do
        read -r sample run path _ <<< "$entry"
        printf '%s.source\t%s/%s_{1,2}.fastq.gz\n' "$sample" "$ENA" "$path"
        printf '%s.raw_pairs\t%s\n'   "$sample" "$(count_pairs "$WORK/raw/${run}_1.fastq.gz")"
        printf '%s.smoke_pairs\t%s\n' "$sample" "$(count_pairs "$WORK/smoke/${sample}_R1.fastq.gz")"
        printf '%s.chr22_chrM_pairs_in_smoke\t%s\n' "$sample" "$(wc -l < "$WORK/align/$sample.hg38.names")"
        printf '%s.ecoli_pairs_in_smoke\t%s\n'      "$sample" "$(wc -l < "$WORK/align/$sample.ecoli.names")"
        printf '%s.ci_pairs\t%s\n'    "$sample" "$(count_pairs "$OUT/fastq/${sample}_R1.fastq.gz")"
    done
    printf 'seed\t%s\n' "$SEED"
    printf 'gtf\t%s\n' "$GTF_URL"
    printf 'ccre\t%s\n' "$CCRE_URL"
    printf 'blacklist\t%s\n' "$BLACKLIST_URL"
    printf 'fasta\t%s/chromosomes/{chr22,chrM}.fa.gz\n' "$UCSC"
    printf 'chrom_sizes\t%s/bigZips/hg38.chrom.sizes\n' "$UCSC"
    printf 'ecoli_fasta\t%s\n' "$ECOLI_URL"
    printf 'selection_index\t%s\n' "$HG38_BT2"
    printf 'bowtie2\t%s\n' "$(bowtie2 --version | head -1 | awk '{ print $NF }')"
    printf 'samtools\t%s\n' "$(samtools --version | head -1 | awk '{ print $2 }')"
    printf 'seqtk\t%s\n' "$(seqtk 2>&1 | awk '/Version/ { print $2 }')"
} > "$OUT/MANIFEST.tsv"

log "done"
column -t -s $'\t' "$OUT/MANIFEST.tsv"
