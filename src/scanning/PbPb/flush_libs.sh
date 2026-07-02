#!/usr/bin/env bash
# Removes stale ACLiC-compiled ROOT libraries from the scan directory.
# Run this before any new scan to force a clean recompile of all macros.
#
# Usage (from anywhere):
#   bash /path/to/flush_libs.sh
# Or, if executable:
#   ./flush_libs.sh

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "Flushing ACLiC libraries in: $DIR"
echo ""

# Patterns produced by ROOT ACLiC (gROOT->LoadMacro("Foo.C+"))
#   Foo_C.so           — compiled shared library
#   Foo_C.d            — dependency file
#   Foo_C_ACLiC_*      — dict headers, .cxx, .o, .pcm, linkdef, map files
#   Foo_C_ACLiC_dict.cxx_tmp_* — leftover temp files from interrupted builds
patterns=(
    "*_C.so"
    "*_C.d"
    "*_C_ACLiC_*"
)

count=0
for pat in "${patterns[@]}"; do
    for f in $pat; do
        if [[ -e "$f" ]]; then
            echo "  rm  $f"
            rm -f "$f"
            (( count++ )) || true
        fi
    done
done

echo ""
if [[ $count -eq 0 ]]; then
    echo "Nothing to remove — directory is already clean."
else
    echo "Removed $count file(s). Next scan will recompile from source."
fi
