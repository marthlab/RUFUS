#!/bin/bash
# Build the small, deterministic fixtures the RUFUS functional suite runs against.
#
# Produces three configurations from one tiny slice of real GRCh38:
#   trio/        child + mother + father  -> a true DE NOVO (child-only) over shared germline
#   somatic/     tumor + normal           -> somatic variants over shared germline
#   specificity/ subjectA + subjectB      -> SAME genome, different seeds; RUFUS must call NOTHING
#
# Everything is seeded, and wgsim runs with -r 0 -R 0 so the ONLY variants present are the ones
# we plant. Re-running reproduces byte-identical fixtures.
#
# Usage:  bash make_fixtures.sh [/path/to/GRCh38.fa]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIX="$HERE/fixtures"

# ---------------- config ----------------
SRC_REF=${1:-/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/resources/reg_test_files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa}
SRC_REGION=chr20:35000000-35199999   # 200 kb of q-arm chr20: N-free, real GC/repeat structure, tiny to run
CONTIG=chr20                          # MUST stay an autosome name -- RUFUS post-filters with VilterAutosomeOnly
LEN=200000
DEPTH=30                              # per-genome coverage; het variants get ~15 reads (comfortably over -m 5)
RL=150
ERR=0.001                             # wgsim base error rate
# ----------------------------------------

# ---- tool resolution (PATH first, then CHPC installdirs) ----
resolve() {  # resolve <name> <glob>...
  local n=$1; shift
  local p; p=$(command -v "$n" 2>/dev/null) && { echo "$p"; return; }
  for g in "$@"; do for c in $g; do [ -x "$c" ] && { echo "$c"; return; }; done; done
  echo "ERROR: cannot find '$n' (set ${n^^}= to override)" >&2; exit 1
}
SAMTOOLS=${SAMTOOLS:-$(resolve samtools '/uufs/chpc.utah.edu/sys/installdir/samtools/*/bin/samtools')}
BCFTOOLS=${BCFTOOLS:-$(resolve bcftools '/uufs/chpc.utah.edu/sys/installdir/bcftools/*/bin/bcftools')}
WGSIM=${WGSIM:-$(resolve wgsim '/uufs/chpc.utah.edu/sys/installdir/samtools/*/bin/wgsim')}
BWA=${BWA:-$(resolve bwa '/uufs/chpc.utah.edu/sys/installdir/bwa/*/bin/bwa' '/uufs/chpc.utah.edu/sys/installdir/bwa/*/bwa')}
BGZIP=${BGZIP:-$(resolve bgzip '/uufs/chpc.utah.edu/sys/installdir/*/bin/bgzip')}
echo "tools: samtools=$SAMTOOLS bcftools=$BCFTOOLS wgsim=$WGSIM bwa=$BWA"

rm -rf "$FIX"; mkdir -p "$FIX"/{ref,designed,trio,somatic,specificity,fastq}

# ============ 1. tiny reference ============
echo "== extracting $SRC_REGION -> tiny reference =="
[ -f "$SRC_REF" ] || { echo "ERROR: source reference not found: $SRC_REF"; exit 1; }
TINY=$FIX/ref/tiny.fa
{ echo ">$CONTIG"; $SAMTOOLS faidx "$SRC_REF" "$SRC_REGION" | tail -n +2; } > "$TINY"
$SAMTOOLS faidx "$TINY"

# refuse an N-heavy slice -- it would silently starve k-mer counting
NFRAC=$(grep -v '^>' "$TINY" | tr -cd 'Nn' | wc -c)
echo "   N bases in slice: $NFRAC / $LEN"
[ "$NFRAC" -lt $((LEN / 20)) ] || { echo "ERROR: slice is >5% N -- pick a different SRC_REGION"; exit 1; }

# ============ 2. designed variants ============
# Built programmatically because REF alleles must match the extracted sequence exactly.
refbase() { $SAMTOOLS faidx "$TINY" "$CONTIG:$1-$2" | tail -n +2 | tr -d '\n' | tr 'acgtn' 'ACGTN'; }
flip()    { case "$1" in A) echo C;; C) echo A;; G) echo T;; T) echo G;; *) echo A;; esac; }

vcf_header() {
  printf '##fileformat=VCFv4.2\n'
  printf '##INFO=<ID=CLASS,Number=1,Type=String,Description="Planted variant class">\n'
  printf '##contig=<ID=%s,length=%d>\n' "$CONTIG" "$LEN"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
}
# emit_snv POS | emit_del POS LEN | emit_ins POS LEN | emit_mnv POS LEN
emit_snv() { local p=$1 r a; r=$(refbase "$p" "$p"); a=$(flip "$r"); printf '%s\t%d\t.\t%s\t%s\t.\t.\tCLASS=SNV\n' "$CONTIG" "$p" "$r" "$a"; }
emit_del() { local p=$1 l=$2 r; r=$(refbase "$p" $((p+l))); printf '%s\t%d\t.\t%s\t%s\t.\t.\tCLASS=DEL%d\n' "$CONTIG" "$p" "$r" "${r:0:1}" "$l"; }
emit_ins() { local p=$1 l=$2 r ins; r=$(refbase "$p" "$p"); ins=$(printf 'ACGTACGTAC' | cut -c1-"$l"); printf '%s\t%d\t.\t%s\t%s%s\t.\t.\tCLASS=INS%d\n' "$CONTIG" "$p" "$r" "$r" "$ins" "$l"; }
emit_mnv() { local p=$1 l=$2 r a i c; r=$(refbase "$p" $((p+l-1))); a=""; for ((i=0;i<${#r};i++)); do c=$(flip "${r:$i:1}"); a="$a$c"; done
             printf '%s\t%d\t.\t%s\t%s\t.\t.\tCLASS=MNV%d\n' "$CONTIG" "$p" "$r" "$a" "$l"; }

echo "== designing variant sets =="
{ vcf_header; emit_snv 20000;  emit_del 40000 10; }                                              > "$FIX/designed/germline.vcf"
{ vcf_header; emit_snv 80000;  emit_ins 90000 5;  emit_mnv 100000 3; }                           > "$FIX/designed/denovo.vcf"
{ vcf_header; emit_snv 120000; emit_del 130000 10; emit_ins 140000 5; emit_mnv 150000 3; emit_del 170000 1000; } > "$FIX/designed/somatic.vcf"

for v in germline denovo somatic; do
  $BGZIP -f -c "$FIX/designed/$v.vcf" > "$FIX/designed/$v.vcf.gz"
  $BCFTOOLS index -f -t "$FIX/designed/$v.vcf.gz"
done
# combined sets used to build genomes (kept on separate lines so `set -e` actually fires on failure)
$BCFTOOLS concat -a "$FIX/designed/germline.vcf.gz" "$FIX/designed/denovo.vcf.gz"  -Oz -o "$FIX/designed/child_all.vcf.gz"
$BCFTOOLS index -f -t "$FIX/designed/child_all.vcf.gz"
$BCFTOOLS concat -a "$FIX/designed/germline.vcf.gz" "$FIX/designed/somatic.vcf.gz" -Oz -o "$FIX/designed/tumor_all.vcf.gz"
$BCFTOOLS index -f -t "$FIX/designed/tumor_all.vcf.gz"

# ============ 3. genomes ============
echo "== building genomes via bcftools consensus =="
$BCFTOOLS consensus -f "$TINY" "$FIX/designed/germline.vcf.gz"  > "$FIX/ref/parent.fa"    # germline background
$BCFTOOLS consensus -f "$TINY" "$FIX/designed/child_all.vcf.gz" > "$FIX/ref/child.fa"     # germline + de novo
$BCFTOOLS consensus -f "$TINY" "$FIX/designed/tumor_all.vcf.gz" > "$FIX/ref/tumor.fa"     # germline + somatic

# ============ 4. reads ============
# npairs for full coverage of LEN at DEPTH with paired RL reads
NFULL=$(( LEN * DEPTH / (2 * RL) ))
NHALF=$(( NFULL / 2 ))
sim() {  # sim <genome.fa> <npairs> <seed> <outprefix>
  $WGSIM -N "$2" -1 $RL -2 $RL -r 0 -R 0 -X 0 -e $ERR -S "$3" "$1" "$4.1.fq" "$4.2.fq" >/dev/null 2>&1
}
echo "== simulating reads (full=$NFULL pairs, half=$NHALF) =="
# trio: parents are pure parent.fa; child is 50/50 parent+child -> de novo is HET, germline is HOM
sim "$FIX/ref/parent.fa" "$NFULL" 11 "$FIX/fastq/mother"
sim "$FIX/ref/parent.fa" "$NFULL" 22 "$FIX/fastq/father"
sim "$FIX/ref/parent.fa" "$NHALF" 31 "$FIX/fastq/child_a"; sim "$FIX/ref/child.fa" "$NHALF" 32 "$FIX/fastq/child_b"
cat "$FIX/fastq/child_a.1.fq" "$FIX/fastq/child_b.1.fq" > "$FIX/fastq/child.1.fq"
cat "$FIX/fastq/child_a.2.fq" "$FIX/fastq/child_b.2.fq" > "$FIX/fastq/child.2.fq"
# somatic: normal is pure germline background; tumor is 50/50 -> somatic HET at VAF~0.5
sim "$FIX/ref/parent.fa" "$NFULL" 41 "$FIX/fastq/normal"
sim "$FIX/ref/parent.fa" "$NHALF" 51 "$FIX/fastq/tumor_a"; sim "$FIX/ref/tumor.fa" "$NHALF" 52 "$FIX/fastq/tumor_b"
cat "$FIX/fastq/tumor_a.1.fq" "$FIX/fastq/tumor_b.1.fq" > "$FIX/fastq/tumor.1.fq"
cat "$FIX/fastq/tumor_a.2.fq" "$FIX/fastq/tumor_b.2.fq" > "$FIX/fastq/tumor.2.fq"
# specificity: SAME genome, different seeds -> expected call set is EMPTY
sim "$FIX/ref/parent.fa" "$NFULL" 61 "$FIX/fastq/subjectA"
sim "$FIX/ref/parent.fa" "$NFULL" 62 "$FIX/fastq/subjectB"
rm -f "$FIX"/fastq/{child_a,child_b,tumor_a,tumor_b}.[12].fq

# ============ 5. align -> BAM + CRAM ============
echo "== indexing tiny reference and aligning =="
$BWA index "$TINY" 2>/dev/null
align() {  # align <name> <outdir>
  local n=$1 d=$2
  $BWA mem -t 4 -R "@RG\tID:$n\tSM:$n\tPL:ILLUMINA\tLB:$n" "$TINY" \
      "$FIX/fastq/$n.1.fq" "$FIX/fastq/$n.2.fq" 2>/dev/null \
    | $SAMTOOLS sort -o "$d/$n.bam" -
  $SAMTOOLS index "$d/$n.bam"
  $SAMTOOLS view -C -T "$TINY" -o "$d/$n.cram" "$d/$n.bam"
  $SAMTOOLS index "$d/$n.cram"
}
for n in child mother father;  do align "$n" "$FIX/trio";        done
for n in tumor normal;         do align "$n" "$FIX/somatic";     done
for n in subjectA subjectB;    do align "$n" "$FIX/specificity"; done

# FASTQs are the bulk of the fixture size; gzip so the set stays committable (~92M -> ~22M).
# The FASTQ-input test can zcat these if RUFUS.Filter wants plain text.
echo "== compressing fastq =="
gzip -f "$FIX"/fastq/*.fq

# ============ 6. expectations + manifest ============
cp "$FIX/designed/denovo.vcf"  "$FIX/trio/expected.vcf"       # child-only; germline must NOT appear
cp "$FIX/designed/somatic.vcf" "$FIX/somatic/expected.vcf"    # tumor-only; germline must NOT appear
{ vcf_header; } > "$FIX/specificity/expected.vcf"             # deliberately EMPTY

{
  echo "# RUFUS functional fixtures -- generated by make_fixtures.sh"
  echo "source_ref=$SRC_REF"
  echo "source_region=$SRC_REGION  (renamed to $CONTIG, length $LEN)"
  echo "depth=${DEPTH}x  read_len=$RL  err=$ERR  wgsim -r 0 -R 0 (only planted variants present)"
  echo
  echo "trio/        child,mother,father   expected.vcf = de novo set (3 variants), germline must NOT be called"
  echo "somatic/     tumor,normal          expected.vcf = somatic set (5 variants), germline must NOT be called"
  echo "specificity/ subjectA,subjectB     expected.vcf = EMPTY -- same genome, different seeds"
  echo
  echo "planted variant classes: SNV, INS5, DEL10, MNV3, DEL1000"
} > "$FIX/MANIFEST.txt"

echo
echo "== done =="
du -sh "$FIX"
cat "$FIX/MANIFEST.txt"
