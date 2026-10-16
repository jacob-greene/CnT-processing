#!/usr/bin/env bash
# Sample sheet validator tests. Needs only Nextflow (no containers, no jobs).
#
#   bash tests/test_samplesheet.sh
#
# Each sheet in tests/samplesheets/ uses the placeholder DATA/ for its FASTQ
# folder. The script creates empty dummy FASTQs in a temporary folder, writes a
# copy of each sheet that points there, and runs the pipeline with
# --validate_only. good*.csv must pass; bad*.csv must fail with the expected message.
set -uo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$tmp/data"
for f in a_L1_R1 a_L1_R2 a_L2_R1 a_L2_R2 b_R1 b_R2 c_R1 c_R2; do
    : | gzip > "$tmp/data/$f.fastq.gz"
done

# sheet -> text that must appear in the error message
declare -A expect=(
    [bad_control_missing.csv]="lists control 'IgG_repX', which is not a sample in this sheet"
    [bad_control_not_control.csv]="does not have is_control=true"
    [bad_control_has_control.csv]="is a control (is_control=true) but itself lists control"
    [bad_conflicting_metadata.csv]="has pool='pool2', but line 2 has pool='pool1'"
    [bad_missing_fastq.csv]="FASTQ not found"
    [bad_r1_r2_mismatch.csv]="R1/R2 count mismatch"
    [bad_is_control_value.csv]="is_control must be 'true' or 'false'"
)

pass=0; fail=0
for sheet in "$repo"/tests/samplesheets/*.csv; do
    name=$(basename "$sheet")
    sed "s#DATA/#$tmp/data/#g" "$sheet" > "$tmp/$name"
    log="$tmp/$name.log"
    (cd "$tmp" && nextflow -q run "$repo/main.nf" --validate_only --genome hg38 \
        --samplesheet "$tmp/$name" > "$log" 2>&1)
    status=$?
    if [[ $name == good* ]]; then
        if [[ $status -eq 0 ]]; then
            echo "PASS  $name (accepted)"; pass=$((pass + 1))
        else
            echo "FAIL  $name was rejected:"; sed 's/^/      /' "$log"; fail=$((fail + 1))
        fi
    else
        if [[ $status -ne 0 ]] && grep -qF "${expect[$name]}" "$log"; then
            echo "PASS  $name (rejected: ${expect[$name]})"; pass=$((pass + 1))
        else
            echo "FAIL  $name (exit $status, expected message: ${expect[$name]})"
            sed 's/^/      /' "$log"; fail=$((fail + 1))
        fi
    fi
done

echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
