#!/usr/bin/env bash
# Compare freshly-rendered viz output against the baseline.
# 5 of 6 plots are bit-deterministic (md5). 05_deletions_binom uses geom_jitter,
# so it is compared perceptually with a pixel-difference tolerance.
#
# Usage: data-raw/viz_compare.sh <current_dir> [baseline_dir]
set -u
CUR="${1:-/tmp/rabhit_viz/current}"
BASE="${2:-/tmp/rabhit_viz/baseline}"
# pixel-diff budget for the jittered plot (well below a real structural change)
JITTER_TOL=20000

status=0
for f in 01_haplotype_single 02_hap_heatmap 03_hap_dendro 04_deletion_heatmap 06_deletions_vpooled; do
  if cmp -s "$BASE/$f.png" "$CUR/$f.png"; then
    echo "OK    $f (identical)"
  else
    echo "DIFF  $f (md5 mismatch) <-- investigate"
    status=1
  fi
done

# jittered plot: perceptual compare
ae=$(compare -metric AE "$BASE/05_deletions_binom.png" "$CUR/05_deletions_binom.png" null: 2>&1 | tr -d '\n')
ae_int=${ae%%.*}
if [[ "$ae_int" =~ ^[0-9]+$ ]] && (( ae_int <= JITTER_TOL )); then
  echo "OK    05_deletions_binom (jitter, ${ae} px diff <= ${JITTER_TOL})"
else
  echo "DIFF  05_deletions_binom (${ae} px diff > ${JITTER_TOL}) <-- investigate"
  status=1
fi

exit $status
