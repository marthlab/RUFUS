#!/bin/bash
# Replay RUFUS.interpret standalone against a preserved run directory.
#
# RUFUS.interpret is the last pipeline stage. Every input it takes is a file that
# earlier stages leave behind in the work directory, so it can be re-run on its own
# in minutes instead of re-running the multi-hour assembly pipeline. This script
# reconstructs the invocation from a preserved run dir, runs it in an isolated
# shadow work dir (the source run dir is never written to), and normalizes the
# output so two runs can be diffed.
#
#   ./replay.sh run    <run_dir> [out_dir]        replay once, print where output landed
#   ./replay.sh freeze <run_dir> <baseline_dir>   replay and store the result as a baseline
#   ./replay.sh check  <run_dir> <baseline_dir>   replay and diff against a stored baseline
#   ./replay.sh twice  <run_dir>                  replay twice and diff, to prove determinism
#
# <run_dir> is the RUFUS work directory of a completed run, e.g.
#   resources/reg_test_files/runs/chr20_m5/rufus_chr20
#
# Exit status: 0 = match / success, 1 = usage or setup error, 2 = baseline mismatch.

set -euo pipefail

RDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INTERPRET="${RUFUS_INTERPRET:-$RDIR/bin/RUFUS.interpret}"
ADDSA="$RDIR/scripts/AddSAtoReadSame.pl"
SAMTOOLS="${SAMTOOLS:-samtools}"

# The reference is not stored in the run dir; override if yours lives elsewhere.
# reg_test_files lives one level above the repo, alongside it under RUFUS/.
REF="${RUFUS_REPLAY_REF:-$RDIR/../resources/reg_test_files/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa}"

die() { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[replay] $*" >&2; }

# Exactly one match for a glob, or die. Keeps a silently-missing or ambiguous
# input from turning into a confusing interpret failure later.
one() {
    local desc="$1"; shift
    local -a hits=()
    local f
    for f in "$@"; do [ -e "$f" ] && hits+=("$f"); done
    [ ${#hits[@]} -eq 0 ] && die "no $desc found (looked for: $*)"
    [ ${#hits[@]} -gt 1 ] && die "$desc is ambiguous, found ${#hits[@]}: ${hits[*]}"
    printf '%s' "${hits[0]}"
}

# ---------------------------------------------------------------------------
# Build an isolated work dir of symlinks mirroring the original layout.
# interpret resolves -hf/-mod/-o relative to CWD and writes outputs (plus six
# auxiliary streams into Intermediates/) under $WORK_DIR, so it needs a writable
# tree. Symlinking keeps the reg-test data pristine -- it is the oracle.
# ---------------------------------------------------------------------------
build_shadow() {
    local run_dir="$1" work="$2"
    rm -rf "$work"
    mkdir -p "$work/Intermediates"

    local f
    for f in "$run_dir"/*; do
        [ -f "$f" ] && ln -sf "$f" "$work/$(basename "$f")"
    done
    for f in "$run_dir"/Intermediates/*; do
        [ -f "$f" ] && ln -sf "$f" "$work/Intermediates/$(basename "$f")"
    done
}

# ---------------------------------------------------------------------------
# Reconstruct the argument list. Mirrors scripts/Overlap.shorter.sh:373.
# ---------------------------------------------------------------------------
run_interpret() {
    local run_dir="$1" out_dir="$2"
    local work="$out_dir/work"

    [ -d "$run_dir" ] || die "run dir not found: $run_dir"
    # Absolute, or the symlinks in the shadow dir resolve relative to the wrong place.
    run_dir="$(cd "$run_dir" && pwd)"
    [ -x "$INTERPRET" ] || die "interpret binary not executable: $INTERPRET"
    [ -r "$REF" ] || die "reference not readable: $REF (set RUFUS_REPLAY_REF)"
    command -v "$SAMTOOLS" >/dev/null || die "samtools not on PATH (set SAMTOOLS)"

    mkdir -p "$out_dir"
    build_shadow "$run_dir" "$work"

    local bam hashlist sample sampleref exclude mob moddist
    bam=$(one       "contigs bam"   "$work"/*.overlap.hashcount.fastq.bam)
    hashlist=$(one  "HashList"      "$work"/*.HashList)
    sample=$(one    "-s sample"     "$work"/Intermediates/*.overlap.asembly.hash.fastq.sample)
    sampleref=$(one "-sR ref sample" "$work"/Intermediates/*.overlap.asembly.hash.fastq.Ref.sample)
    exclude=$(one   "-e RepRefHash" "$work"/Intermediates/*.ref.RepRefHash)
    mob=$(one       "-mob MOB.sam"  "$work"/Intermediates/*.MOB.sam)

    # The model file is frequently absent (see docs/RUFUS.interpret.audit.md 2.1).
    # Pass it regardless -- a faithful replay must reproduce that, and the log
    # scrape below turns its absence into a visible, diffable signal.
    # Derive it from the subject stub the same way Overlap.shorter.sh does (split the
    # NameStub on ".V2"); globbing would pick up a control's histogram instead.
    local subj_stub; subj_stub="$(basename "$bam" .V2.overlap.hashcount.fastq.bam)"
    moddist="$work/$subj_stub.Jhash.histo.7.7.dist"

    # Controls: -c and -cR must be paired in matching order. Sorting both lists by
    # the control name they end with keeps them aligned (interpret pairs by ordinal
    # position only, with no consistency check -- audit Tier 2, B14).
    local -a ctrl_args=()
    local c cr
    while IFS= read -r c; do
        [ -n "$c" ] || continue
        cr="$(dirname "$c")/$(basename "$c" | sed 's/^ctrlhash\.overlap\.asembly\.hash\.fastq\./ctrlhash.overlap.asembly.hash.fastq.Ref./')"
        [ -e "$cr" ] || die "control $c has no matching .Ref. file at $cr"
        ctrl_args+=(-c "$c" -cR "$cr")
    done < <(find "$work/Intermediates" -maxdepth 1 -name 'ctrlhash.overlap.asembly.hash.fastq.*' \
                  ! -name '*.fastq.Ref.*' | sort)

    [ ${#ctrl_args[@]} -eq 0 ] && note "warning: no control hashes found; running with no -c"

    local stub; stub="$(basename "$bam")"

    note "run dir : $run_dir"
    note "work    : $work"
    note "controls: $(( ${#ctrl_args[@]} / 4 ))"
    [ -e "$moddist" ] && note "model   : present" || note "model   : ABSENT ($(basename "$moddist"))"

    # -rp points at the repo root so interpret can find resources/vcf_header.txt.
    # Relative -hf/-mod/-o require CWD == WORK_DIR.
    local rc=0
    (
        cd "$work"
        WORK_DIR="$work" "$SAMTOOLS" view -h "$bam" \
            | perl "$ADDSA" \
            | grep -v chrUn \
            | WORK_DIR="$work" "$INTERPRET" \
                -mob "$mob" \
                -mod "$(basename "$moddist")" \
                -mQ 10 \
                -r "$REF" \
                -hf "$(basename "$hashlist")" \
                -o "$stub" \
                -m 1000 \
                "${ctrl_args[@]}" \
                -sR "$sampleref" \
                -s "$sample" \
                -e "$exclude" \
                -rp "$RDIR" \
                -plct 7
    ) > "$out_dir/interpret.stdout.log" 2>"$out_dir/interpret.stderr.log" || rc=$?

    echo "$rc" > "$out_dir/exit_code"
    note "exit code: $rc"

    local vcf="$work/$stub.vcf"
    [ -f "$vcf" ] || die "interpret produced no VCF at $vcf (see $out_dir/interpret.stdout.log)"

    normalize_vcf "$vcf" > "$out_dir/calls.vcf"
    scrape_log "$out_dir/interpret.stdout.log" "$out_dir/interpret.stderr.log" > "$out_dir/signals.txt"
    summarize "$out_dir/calls.vcf" > "$out_dir/summary.txt"

    note "wrote $out_dir/{calls.vcf,signals.txt,summary.txt,exit_code}"
}

# ---------------------------------------------------------------------------
# Strip the two header lines that change on every run for reasons unrelated to
# behaviour: the epoch timestamp and the command line (which embeds absolute
# paths). Everything else is kept, in order -- record order is meaningful.
# ---------------------------------------------------------------------------
normalize_vcf() {
    grep -v '^##fileDate=' "$1" | grep -v '^##RUFUSCommandLine='
}

# ---------------------------------------------------------------------------
# interpret prints megabytes of debug to stdout and returns 0 even on fatal
# errors (audit Tier 3), so neither the full log nor the exit code is a usable
# regression signal on its own. Pull out the lines that indicate a load failure
# or an internal complaint, dedupe with counts, and diff that instead.
# ---------------------------------------------------------------------------
scrape_log() {
    cat "$@" 2>/dev/null \
        | grep -aiE 'error|cannot|could not|no model file|this is going to break|well shit|WTF|out of sync|no base works|ERROR IN RevComp|not found' \
        | sed 's/[0-9]\{4,\}/<N>/g' \
        | sort | uniq -c | sort -rn
}

summarize() {
    local vcf="$1"
    echo "records:   $(grep -vc '^#' "$vcf" || true)"
    echo
    echo "FILTER:"
    grep -v '^#' "$vcf" | cut -f7 | sort | uniq -c | sort -rn | sed 's/^/  /'
    echo
    echo "GT:"
    grep -v '^#' "$vcf" | awk -F'\t' 'NF>=10{split($10,a,":"); print a[1]}' | sort | uniq -c | sort -rn | sed 's/^/  /'
    echo
    echo "QUAL:"
    grep -v '^#' "$vcf" | awk -F'\t' '
        $6=="."{d++; next}
        {n++; s+=$6; if($6>mx)mx=$6; if($6>100)over++}
        END{printf "  n=%d  mean=%.2f  max=%.4f  >100=%d  missing=%d\n", n, (n?s/n:0), mx, over+0, d+0}'
    echo
    echo "columns (should be constant):"
    grep -v '^#' "$vcf" | awk -F'\t' '{print NF}' | sort | uniq -c | sed 's/^/  /'
}

compare() {
    local out_dir="$1" baseline="$2" status=0

    for f in calls.vcf signals.txt summary.txt exit_code; do
        [ -f "$baseline/$f" ] || die "baseline is missing $f -- re-run 'freeze'"
        [ -f "$out_dir/$f" ]  || die "current run produced no $f -- the replay itself failed"
    done

    echo "=============================================================="
    echo " replay check"
    echo "   baseline: $baseline"
    echo "   current : $out_dir"
    echo "=============================================================="
    echo

    if diff -q "$baseline/calls.vcf" "$out_dir/calls.vcf" >/dev/null; then
        echo "VCF       : identical"
    else
        status=2
        local added removed
        added=$(comm -13 <(grep -v '^#' "$baseline/calls.vcf" | sort) \
                         <(grep -v '^#' "$out_dir/calls.vcf" | sort) | wc -l)
        removed=$(comm -23 <(grep -v '^#' "$baseline/calls.vcf" | sort) \
                           <(grep -v '^#' "$out_dir/calls.vcf" | sort) | wc -l)
        echo "VCF       : CHANGED  (+$added records, -$removed records)"
        echo "            full diff: diff $baseline/calls.vcf $out_dir/calls.vcf"
    fi

    if diff -q "$baseline/signals.txt" "$out_dir/signals.txt" >/dev/null; then
        echo "log signals: identical"
    else
        status=2
        echo "log signals: CHANGED"
        diff "$baseline/signals.txt" "$out_dir/signals.txt" | sed 's/^/            /' || true
    fi

    if [ "$(cat "$baseline/exit_code")" = "$(cat "$out_dir/exit_code")" ]; then
        echo "exit code : unchanged ($(cat "$out_dir/exit_code"))"
    else
        status=2
        echo "exit code : CHANGED  $(cat "$baseline/exit_code") -> $(cat "$out_dir/exit_code")"
    fi

    echo
    echo "--- summary diff (baseline -> current) ---"
    diff "$baseline/summary.txt" "$out_dir/summary.txt" || true

    return $status
}

usage() { sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1; }

cmd="${1:-}"; shift || usage
case "$cmd" in
    run)
        [ $# -ge 1 ] || usage
        run_interpret "$1" "${2:-$(mktemp -d -t rufus-replay-XXXXXX)}"
        ;;
    freeze)
        [ $# -eq 2 ] || usage
        tmp=$(mktemp -d -t rufus-replay-XXXXXX)
        run_interpret "$1" "$tmp"
        mkdir -p "$2"
        cp "$tmp"/{calls.vcf,signals.txt,summary.txt,exit_code} "$2"/
        note "baseline frozen in $2"
        note "commit it: git add $2 && git commit -m 'Freeze interpret replay baseline'"
        ;;
    check)
        [ $# -eq 2 ] || usage
        tmp=$(mktemp -d -t rufus-replay-XXXXXX)
        run_interpret "$1" "$tmp"
        compare "$tmp" "$2"
        ;;
    twice)
        [ $# -eq 1 ] || usage
        a=$(mktemp -d -t rufus-replay-a-XXXXXX); b=$(mktemp -d -t rufus-replay-b-XXXXXX)
        run_interpret "$1" "$a"
        run_interpret "$1" "$b"
        if diff -q "$a/calls.vcf" "$b/calls.vcf" >/dev/null && \
           diff -q "$a/signals.txt" "$b/signals.txt" >/dev/null; then
            echo "DETERMINISTIC: two runs produced identical output."
        else
            echo "NON-DETERMINISTIC: runs differ. A baseline cannot be trusted until this is fixed."
            diff "$a/calls.vcf" "$b/calls.vcf" | head -40 || true
            exit 2
        fi
        ;;
    *) usage ;;
esac
