#!/bin/sh

# Layout-drift guard for the prepfold / prepfold_multi shared-struct trick.
#
# prepfold_multi reuses prepfold_pipeline.c, so its parsed options must populate
# the SAME `struct s_Cmdline` memory layout that prepfold uses.  This is achieved
# by making include/prepfold_multi_cmd.h a copy of include/prepfold_cmd.h with the
# extra `candfile` field appended AFTER `full_cmd_line` -- every field the shared
# pipeline touches therefore keeps prepfold's offset.
#
# This script enforces that invariant structurally: it extracts the body of
# `typedef struct s_Cmdline { ... } Cmdline;` from both headers and asserts that
# prepfold_cmd.h's body is a line-for-line PREFIX of prepfold_multi_cmd.h's body
# (the only difference being the extra trailing -candfile block).  Any reorder,
# rename, type change, or removal of a shared field breaks the prefix and fails
# the check -- catching silent ABI drift before it corrupts shared-struct reads.
#
# Exits 0 on PASS, 1 on drift, 2 on a usage/parse problem.  BASE and MULTI may be
# overridden (env vars) to point at alternate headers -- used by the self-test.

_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INC="${_SCRIPT_DIR}/../../include"
BASE="${BASE:-${INC}/prepfold_cmd.h}"
MULTI="${MULTI:-${INC}/prepfold_multi_cmd.h}"

# Print the struct body (lines strictly between the opening brace line and the
# closing `} Cmdline;` line) of the given header.
extract_body() {
    awk '
        /typedef[ \t]+struct[ \t]+s_Cmdline[ \t]*\{/ { grab = 1; next }
        grab && /^[ \t]*\}[ \t]*Cmdline[ \t]*;/      { grab = 0 }
        grab                                         { print }
    ' "$1"
}

for f in "$BASE" "$MULTI"; do
    if [ ! -f "$f" ]; then
        echo "FAIL: cmd-layout: header not found: $f"
        exit 2
    fi
done

tmpdir="$(mktemp -d)" || { echo "FAIL: cmd-layout: mktemp failed"; exit 2; }
trap 'rm -rf "$tmpdir"' EXIT INT TERM

extract_body "$BASE"  > "$tmpdir/base"
extract_body "$MULTI" > "$tmpdir/multi"

nbase=$(wc -l < "$tmpdir/base")
nmulti=$(wc -l < "$tmpdir/multi")

if [ "$nbase" -eq 0 ] || [ "$nmulti" -eq 0 ]; then
    echo "FAIL: cmd-layout: could not extract s_Cmdline body (header format drift?)"
    exit 2
fi

# prepfold_multi must add fields, never remove: its body has to be longer.
if [ "$nmulti" -le "$nbase" ]; then
    echo "FAIL: cmd-layout: prepfold_multi_cmd.h body ($nmulti lines) is not longer"
    echo "      than prepfold_cmd.h body ($nbase lines) -- the appended -candfile"
    echo "      block is missing or shared fields were removed."
    exit 1
fi

head -n "$nbase" "$tmpdir/multi" > "$tmpdir/multi_prefix"

if diff "$tmpdir/base" "$tmpdir/multi_prefix" >/dev/null 2>&1; then
    extra=$((nmulti - nbase))
    echo "PASS: cmd-layout (prepfold_cmd.h is a prefix of prepfold_multi_cmd.h; +${extra} trailing lines)"
    exit 0
fi

echo "FAIL: cmd-layout: prepfold_cmd.h's s_Cmdline is NOT a prefix of"
echo "      prepfold_multi_cmd.h's -- shared prepfold_pipeline.c field offsets may"
echo "      have drifted.  Diff (prepfold_cmd.h '<' vs prepfold_multi_cmd.h prefix '>'):"
diff "$tmpdir/base" "$tmpdir/multi_prefix" | sed 's/^/        /'
exit 1
