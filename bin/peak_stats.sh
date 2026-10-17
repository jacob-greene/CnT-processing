#!/usr/bin/env bash
# Peak stats and FRiP for one target: one header line + one row of peaks_frip.tsv.
#
# Usage: peak_stats.sh <id> <control id|NA> <seacr mode> <fragments.bed.gz> <peaks.bed> \
#                      <control fragments.bed.gz | ""> <top0.01 peaks.bed | "">
#
#   peaks          number of peaks
#   median_width, mean_width, bp_in_peaks   from peak end - start
#   frip           target fragments overlapping >= 1 peak / all target fragments
#   control_frip   control fragments overlapping >= 1 of these peaks / all control fragments
# The top01_* columns repeat this for the SEACR 0.01 set when a control exists
# (otherwise the primary set IS the 0.01 set and they are NA).
set -euo pipefail
id=$1; control=$2; mode=$3; frags=$4; peaks=$5; ctrl_frags=${6:-}; top=${7:-}

# width stats: "<n> <median> <mean> <bp>"
widths() {
    awk '{ print $3 - $2 }' "$1" | sort -n | awk '
        { w[NR] = $1; s += $1 }
        END { if (NR == 0) { print "0 NA NA 0"; exit }
              m = NR % 2 ? w[(NR + 1) / 2] : (w[NR / 2] + w[NR / 2 + 1]) / 2
              printf "%d %s %.1f %d\n", NR, m, s / NR, s }'
}
# fragments overlapping >= 1 peak, and all fragments: "<in> <total>"
in_peaks() {
    local total n
    total=$(zcat "$1" | wc -l)
    if [[ -s "$2" ]]; then
        n=$(bedtools intersect -u -a "$1" -b "$2" | wc -l)
    else
        n=0
    fi
    echo "$n $total"
}
frac() { awk -v a="$1" -v b="$2" 'BEGIN { if (b > 0) printf "%.6f", a / b; else print "NA" }'; }

read -r n med mean bp <<< "$(widths "$peaks")"
read -r f_in f_tot <<< "$(in_peaks "$frags" "$peaks")"
c_in=NA; c_tot=NA; c_frip=NA
if [[ -n "$ctrl_frags" ]]; then
    read -r c_in c_tot <<< "$(in_peaks "$ctrl_frags" "$peaks")"
    c_frip=$(frac "$c_in" "$c_tot")
fi
t_n=NA; t_med=NA; t_mean=NA; t_bp=NA; t_frip=NA; t_cfrip=NA
if [[ -n "$top" ]]; then
    read -r t_n t_med t_mean t_bp <<< "$(widths "$top")"
    read -r t_in _ <<< "$(in_peaks "$frags" "$top")"
    t_frip=$(frac "$t_in" "$f_tot")
    if [[ -n "$ctrl_frags" ]]; then
        read -r tc_in _ <<< "$(in_peaks "$ctrl_frags" "$top")"
        t_cfrip=$(frac "$tc_in" "$c_tot")
    fi
fi

printf 'sample\tcontrol\tseacr_mode\tpeaks\tmedian_width\tmean_width\tbp_in_peaks\tfragments\tfragments_in_peaks\tfrip\tcontrol_fragments\tcontrol_fragments_in_peaks\tcontrol_frip\ttop01_peaks\ttop01_median_width\ttop01_mean_width\ttop01_bp_in_peaks\ttop01_frip\ttop01_control_frip\n'
printf '%s\t' "$id" "$control" "$mode" "$n" "$med" "$mean" "$bp" "$f_tot" "$f_in" "$(frac "$f_in" "$f_tot")" \
    "$c_tot" "$c_in" "$c_frip" "$t_n" "$t_med" "$t_mean" "$t_bp" "$t_frip"
printf '%s\n' "$t_cfrip"
