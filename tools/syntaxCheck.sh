#!/usr/bin/env bash
# Type-check a scan macro locally, including the #ifdef DO_FASTJET block.
#
# WHY. The FastJet block is the highest-risk code in this repo to edit and, until
# FastJet was installed here, could not be compiled outside lxplus at all. Scan
# edits went out verified only by brace-balance and eyeball; a typo cost a full
# round trip of hours. This closes that.
#
# Usage:
#   tools/syntaxCheck.sh                                  # PbPb_pfCandAnalyzer.C
#   tools/syntaxCheck.sh src/scanning/pp/pp_scan.C
#
# Exit 0 = clean. This is -fsyntax-only: it type-checks but does not link or
# run, so it catches typos, wrong argument counts, undeclared or out-of-scope
# symbols and include-order mistakes -- not logic errors and not anything that
# only shows up at link time.
#
# The generated driver adds `using namespace std;` because ROOT's interpreter
# and ACLiC have it in effect, and eventMap.h relies on that (it writes bare
# `vector<Float_t>`). Without it g++ reports dozens of spurious errors in
# eventMap.h that ROOT never sees.

set -u

ROOT_BIN="${ROOT_BIN:-$HOME/Programs/root/bin}"
FASTJET_INC="${FASTJET_INC:-$HOME/Programs/fastjet-3.4.3/include}"

export PATH="$ROOT_BIN:$PATH"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${1:-src/scanning/PbPb/PbPb_pfCandAnalyzer.C}"

if [ ! -f "$REPO/$TARGET" ]; then
  echo "ERROR: no such file: $REPO/$TARGET" >&2
  exit 2
fi
if ! command -v root-config >/dev/null 2>&1; then
  echo "ERROR: root-config not found. Set ROOT_BIN." >&2
  exit 2
fi
if [ ! -d "$FASTJET_INC/fastjet" ]; then
  echo "ERROR: FastJet headers not found at $FASTJET_INC" >&2
  echo "  Install with:" >&2
  echo "    cd ~/Programs && curl -LO https://fastjet.fr/repo/fastjet-3.4.3.tar.gz" >&2
  echo "    tar xzf fastjet-3.4.3.tar.gz && cd fastjet-3.4.3" >&2
  echo "    ./configure --prefix=\$HOME/Programs/fastjet" >&2
  echo "  configure alone is enough -- it generates fastjet/config.h, and" >&2
  echo "  -fsyntax-only needs headers only, not the built library." >&2
  exit 2
fi

DIR="$(dirname "$REPO/$TARGET")"
BASE="$(basename "$TARGET")"
DRV="$DIR/._syntaxCheck_$$.C"

cleanup(){ rm -f "$DRV"; }
trap cleanup EXIT

cat > "$DRV" <<EOF
#include <vector>
#include <string>
#include <iostream>
using namespace std;
#include "$BASE"
EOF

echo "checking $TARGET  (DO_FASTJET enabled)"
cd "$DIR" || exit 2
g++ -fsyntax-only -std=c++17 -DDO_FASTJET \
    $(root-config --cflags) \
    -I"$FASTJET_INC" \
    "$(basename "$DRV")" 2>&1 | grep -v "^In file included\|^ *from "

rc=${PIPESTATUS[0]}
if [ "$rc" -eq 0 ]; then echo "OK - no errors"; else echo "FAILED (exit $rc)"; fi
exit "$rc"
