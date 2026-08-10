#!/bin/bash
# tests/accuracy/compare.sh
#
# Score the chr20 accuracy run against the SMaHT COLO829 SNV truth set and print a PASS/FAIL verdict.
# Paths, image, and threshold come from config.sh. Run AFTER run_chr20.slurm finishes:
#
#     cd tests/accuracy && bash compare.sh          # override threshold: MIN_RECALL=0.90 bash compare.sh
#
# Exit status is 0 only if recall meets the threshold.
#   recall    = fraction of TIER1 truth SNVs RUFUS found  <-- the gate
#   precision = fraction of RUFUS SNVs matching ANY truth tier (lower bound; truth isn't exhaustive)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="$HERE/config.sh"
[ -f "$CFG" ] || { echo "ERROR: missing $CFG — copy config.sh.example to config.sh (see README.md)"; exit 1; }
# shellcheck source=/dev/null
source "$CFG"
MIN_RECALL=${MIN_RECALL:-0.85}

WORK=$RUFUS_TEST_DATA/runs/_compare; mkdir -p "$WORK"
module load apptainer 2>/dev/null || true
bcf() { apptainer --quiet exec --bind "$RUFUS_TEST_DATA" "$RUFUS_SIF" bcftools "$@"; }   # --quiet: no gocryptfs chatter

# Tier1 truth (recall denominator) + all-tier truth (precision denominator), restricted to the region.
TT=$WORK/truth.tier1.${REGION}.vcf.gz
bcf view -r "$REGION" -i 'INFO/RGN_T="Tier1"' "$TRUTH_VCF" -Oz -o "$TT"; bcf index -t "$TT"
TA=$WORK/truth.all.${REGION}.vcf.gz
bcf view -r "$REGION" "$TRUTH_VCF" -Oz -o "$TA"; bcf index -t "$TA"
echo "Tier1 truth SNVs on $REGION: $(bcf view -H "$TT" | wc -l)"

echo
echo "recall    = fraction of TIER1 truth SNVs RUFUS found   <-- the GATE"
echo "precision = fraction of RUFUS SNVs that match ANY truth tier   (lower bound; truth isn't exhaustive)"
echo "PASS threshold: recall >= $MIN_RECALL"
echo
overall_pass=1; nrun=0; npass=0
for d in "$RUFUS_TEST_DATA"/runs/chr20_m*/; do
  M=$(basename "$d" | sed 's/chr20_m//')
  nrun=$((nrun+1))
  CALLS=$(ls "$d"/temp.RUFUS.Final.*.vcf.gz 2>/dev/null | head -1) \
    || { printf '  %-4s  no final VCF found%*s[FAIL]\n' "m$M" 44 ""; overall_pass=0; continue; }

  # normalize + SNV-only RUFUS calls (region via -t: streaming filter; -r would need an index -> bgzip error)
  CN=$WORK/calls.m${M}.snv.vcf.gz
  bcf norm -f "$REF_FASTA" "$CALLS" 2>/dev/null | bcf view -v snps -t "$REGION" -Oz -o "$CN" -
  bcf index -t "$CN"

  # RECALL: RUFUS vs TIER1 truth   (0000 = Tier1 missed / FN,  0002 = Tier1 found / TP)
  ISEC=$WORK/isec.m${M}; rm -rf "$ISEC"; bcf isec -p "$ISEC" "$TT" "$CN"
  FN=$(grep -vc '^#' "$ISEC/0000.vcf"); TP=$(grep -vc '^#' "$ISEC/0002.vcf")
  rden=$((TP+FN)); rec=$(awk -v t=$TP -v d=$rden 'BEGIN{printf d?"%.3f":"NA",t/d}')

  # PRECISION: RUFUS vs ALL-tier truth   (0001 = RUFUS in no truth tier / FP,  0002 = matched / TP)
  ISECA=$WORK/isecall.m${M}; rm -rf "$ISECA"; bcf isec -p "$ISECA" "$TA" "$CN"
  FP=$(grep -vc '^#' "$ISECA/0001.vcf"); TPa=$(grep -vc '^#' "$ISECA/0002.vcf")
  pden=$((TPa+FP)); pre=$(awk -v t=$TPa -v d=$pden 'BEGIN{printf d?"%.3f":"NA",t/d}')

  if awk -v r="$rec" -v t="$MIN_RECALL" 'BEGIN{exit !(r!="NA" && r+0>=t+0)}'; then
    verdict="PASS"; npass=$((npass+1))
  else
    verdict="FAIL"; overall_pass=0
  fi
  printf '  %-4s  recall %s/%s = %-6s   precision %s/%s = %-6s   [%s]\n' \
         "m$M" "$TP" "$rden" "$rec" "$TPa" "$pden" "$pre" "$verdict"
done

echo
if [ "$overall_pass" = 1 ]; then
  echo "==================== VERDICT: PASS  ($npass/$nrun met recall >= $MIN_RECALL) ===================="
else
  echo "==================== VERDICT: FAIL  ($npass/$nrun met recall >= $MIN_RECALL) ===================="
fi
echo "(intermediate isec VCFs under: $WORK)"
exit $([ "$overall_pass" = 1 ] && echo 0 || echo 1)
