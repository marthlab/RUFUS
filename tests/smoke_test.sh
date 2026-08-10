#!/bin/bash
# RUFUS container smoke test.
#
# Verifies the image is internally consistent: every binary RUFUS depends on responds, the
# RUFUS executables and run script are present, and the runtime environment is configured.
# It does NOT call variants -- functional/variant-calling tests are a separate suite (TBD).
# Used both as the final build step in the Dockerfile and as the CI gate before pushing.
#
# Exits non-zero on the first failure.
set -uo pipefail

: "${RUFUS_ROOT:=/opt/RUFUS}"
JELLYFISH="${RUFUS_ROOT}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/bin/jellyfish"

fail=0

# check_responds <label> <command...> : passes if the command runs without "command not found".
# Many of these tools exit non-zero when given no args / --version, which is fine — we only
# care that the binary exists and is executable.
check_responds() {
    local label="$1"; shift
    if "$@" >/dev/null 2>&1 || [ $? -ne 127 ]; then
        echo "  ok: ${label}"
    else
        echo "  FAIL: ${label} (not found / not executable)"
        fail=1
    fi
}

check_file() {
    local label="$1" path="$2"
    if [ -e "$path" ]; then
        echo "  ok: ${label} (${path})"
    else
        echo "  FAIL: ${label} missing (${path})"
        fail=1
    fi
}

echo "== external tools =="
check_responds "samtools"  samtools --version
check_responds "bcftools"  bcftools --version
check_responds "bedtools"  bedtools --version
check_responds "bamtools"  bamtools --version
check_responds "bgzip"     bgzip --version
check_responds "aws"       aws --version
check_responds "parallel"  parallel --version

# `bcftools --version` above succeeds even when plugin loading is broken, so check a plugin
# actually dlopens. RUFUS's VCF post-processing uses fill-from-fasta; if BCFTOOLS_PLUGINS points
# at a mismatched (e.g. host-inherited) build, this fails with `undefined symbol: ...` and every
# run dies at the vcf_processing stage. Guard it here so CI catches it in seconds.
echo "== bcftools plugins =="
if [ -n "${BCFTOOLS_PLUGINS:-}" ]; then
    echo "  ok: BCFTOOLS_PLUGINS set (${BCFTOOLS_PLUGINS})"
else
    echo "  FAIL: BCFTOOLS_PLUGINS not set — host env can hijack plugin loading"
    fail=1
fi
if bcftools +fill-from-fasta --version >/dev/null 2>&1; then
    echo "  ok: fill-from-fasta plugin loads"
else
    echo "  FAIL: fill-from-fasta plugin will not load (check BCFTOOLS_PLUGINS / ABI mismatch)"
    fail=1
fi

echo "== RUFUS binaries =="
check_responds "RUFUS.Filter" RUFUS.Filter
check_responds "ModelDist"    ModelDist
check_file     "modified jellyfish" "$JELLYFISH"

echo "== RUFUS layout / environment =="
check_file "runRufus.sh"  "${RUFUS_ROOT}/runRufus.sh"
check_file "globals.txt"  "${RUFUS_ROOT}/resources/globals.txt"
# Provenance stamp must be present so every image self-identifies (its VALUE may be "unknown"
# for a bare local build; here we only assert the file exists, i.e. the stamp step ran).
check_file "BUILD_INFO"   "${RUFUS_ROOT}/BUILD_INFO"
if [ -n "${RUFUS_ROOT:-}" ]; then
    echo "  ok: RUFUS_ROOT set (${RUFUS_ROOT})"
else
    echo "  FAIL: RUFUS_ROOT not set"
    fail=1
fi

echo
if [ "$fail" -ne 0 ]; then
    echo "SMOKE TEST FAILED"
    exit 1
fi
echo "SMOKE TEST PASSED"
