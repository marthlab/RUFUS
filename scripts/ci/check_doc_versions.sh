#!/bin/bash
# Keep the version strings that must appear literally in the docs in step with
# resources/globals.txt.
#
#   check_doc_versions.sh            verify (exit 1 on drift) -- this is what CI runs
#   check_doc_versions.sh --write    rewrite the managed spots to match globals.txt
#
# Most places that need the version can just derive it (docs/publish_new_sif.md reads
# globals.txt inline, and the workflow extracts it the same way). Only user-facing text
# that is read outside a checkout -- i.e. README.md -- has to carry a literal, so only
# those spots are managed here. Each rule is an explicit anchor plus the expected
# rendering, so this can never false-positive on a deliberate historical reference to an
# older version.
#
# The README's *download instructions* no longer carry a version at all: the Docker Hub
# route pulls :latest (or a version the reader chooses), and the Zenodo route asks the API
# which file to fetch from the concept record. Both are drift-free by construction, so the
# only literal left is the human-readable tagline. Note that a stable Zenodo asset name
# would NOT have been sufficient on its own -- Zenodo 404s on /records/<concept>/files/...,
# because it redirects the record endpoint but not paths beneath it, so the version-specific
# record id is unavoidable in a hand-built URL. Hence the API lookup.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

GLOBALS="resources/globals.txt"
[ -r "$GLOBALS" ] || { echo "ERROR: cannot read $GLOBALS" >&2; exit 1; }

# Same extraction the workflow's release guard uses -- keep these identical.
VERSION=$(grep -E '^RUFUS_VERSION=' "$GLOBALS" | cut -d'"' -f2)
[ -n "$VERSION" ] || { echo "ERROR: RUFUS_VERSION not found or empty in $GLOBALS" >&2; exit 1; }

MODE="${1:-check}"
case "$MODE" in
    ""|check) MODE=check ;;
    --write)  MODE=write ;;
    *) echo "usage: $0 [--write]" >&2; exit 1 ;;
esac

status=0

# check_rule <file> <anchor-regex> <expected-line> <sed-substitution>
#   anchor-regex   selects exactly the line(s) this rule owns
#   expected-line  what that line must look like once correct
#   sed-subst      how --write repairs it
check_rule() {
    local file="$1" anchor="$2" expected="$3" subst="$4"
    local found
    found=$(grep -nE "$anchor" "$file" || true)

    if [ -z "$found" ]; then
        echo "  FAIL $file: no line matching /$anchor/ -- the doc was restructured, update this rule" >&2
        status=1
        return
    fi

    if [ "$(printf '%s\n' "$found" | wc -l)" -ne 1 ]; then
        echo "  FAIL $file: /$anchor/ matched more than one line; the rule is ambiguous" >&2
        printf '%s\n' "$found" | sed 's/^/         /' >&2
        status=1
        return
    fi

    local line="${found#*:}"
    if [ "$line" = "$expected" ]; then
        echo "  ok   $file: $expected"
        return
    fi

    if [ "$MODE" = write ]; then
        sed -i "$subst" "$file"
        echo "  wrote $file: $line -> $expected"
    else
        echo "  FAIL $file" >&2
        echo "         have: $line" >&2
        echo "         want: $expected" >&2
        status=1
    fi
}

echo "RUFUS_VERSION in $GLOBALS = $VERSION"

# README tagline, e.g. "K-mer based variant detection. v1.2.0."
check_rule README.md \
    '^K-mer based variant detection\.' \
    "K-mer based variant detection. ${VERSION}." \
    "s/^K-mer based variant detection\..*/K-mer based variant detection. ${VERSION}./"

# README "pin a specific version" example. The :latest line beside it is deliberately not
# managed -- it must stay literal.
check_rule README.md \
    '^apptainer pull rufus\.sif docker://stefinfection/rufus:v' \
    "apptainer pull rufus.sif docker://stefinfection/rufus:${VERSION}" \
    "s|^apptainer pull rufus\.sif docker://stefinfection/rufus:v.*|apptainer pull rufus.sif docker://stefinfection/rufus:${VERSION}|"

if [ "$status" -ne 0 ]; then
    cat >&2 <<EOF

Documentation version strings are out of step with $GLOBALS.
Fix with:  bash scripts/ci/check_doc_versions.sh --write
EOF
    exit 1
fi

echo "All managed doc version strings match."
