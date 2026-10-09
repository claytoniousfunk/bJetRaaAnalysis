#!/usr/bin/env bash
# Run pp_jetTrkMax_scan.C over every legacy skim group, in parallel, then hadd.
# Run ON LXPLUS from src/scanning/pp/ (2026-10-09). The skims are small (~77k
# events, ~3 MB each), so this is a local loop rather than condor.
#   ./run_pp_jetTrkMax_scan.sh [nParallel=8] [firstN=all]   # firstN: trial on the first N groups
# (The group list is GRPS: bash's GROUPS is a read-only builtin holding the user's
# group ids, and assigning to it silently does nothing -- the first version of this
# script ran a single bogus group because of that.)
set -uo pipefail
NP=${1:-8}; FIRST=${2:-0}
SKIMS=/eos/cms/store/group/phys_heavyions/cbennett/skims/output_skims_pp_HighEGJet_withJetTrackMaxInfo
OUT=$(grep -o '/eos/[^"]*output_pp_jetTrkMax[^"/]*' pp_jetTrkMax_scan.C | head -1)   # the macro's own output dir
[ -n "$OUT" ] || { echo "ERROR: no output dir found in pp_jetTrkMax_scan.C"; exit 1; }
mkdir -p "$OUT" logs_jetTrkMax || { echo "ERROR: cannot create $OUT"; exit 1; }
GRPS=$(ls "$SKIMS" | sed -n 's/^pp_skim_output_\([0-9]*\)\.root$/\1/p' | sort -n)
[ -n "$GRPS" ] || { echo "ERROR: no skims found in $SKIMS"; exit 1; }
if [ "$FIRST" -gt 0 ]; then GRPS=$(echo "$GRPS" | head -n "$FIRST"); fi
NG=$(echo "$GRPS" | wc -l)
echo "skims: $SKIMS"
echo "output: $OUT"
echo "groups: $NG (first $(echo "$GRPS" | head -1), last $(echo "$GRPS" | tail -1)), $NP in parallel; logs in logs_jetTrkMax/g<group>.log"
# interpreted, as the 2025 scan (ACLiC would lack the implicit `using namespace std`
# that eventMap.h relies on)
echo "$GRPS" | xargs -P "$NP" -I{} sh -c "root -l -b -q 'pp_jetTrkMax_scan.C({})' > logs_jetTrkMax/g{}.log 2>&1 && echo 'group {} ok' || echo 'group {} FAILED (see logs_jetTrkMax/g{}.log)'"
NOUT=$(ls "$OUT"/pp_scan_output_*.root 2>/dev/null | wc -l)
echo "outputs: $NOUT of $NG"
[ "$NOUT" -gt 0 ] || { echo "ERROR: no outputs; first log:"; head -30 logs_jetTrkMax/g$(echo "$GRPS" | head -1).log; exit 1; }
MERGED="$(dirname "$OUT")/$(basename "$OUT").root"
hadd -f "$MERGED" "$OUT"/pp_scan_output_*.root > logs_jetTrkMax/hadd.log 2>&1 && echo "merged: $MERGED" || echo "hadd FAILED, see logs_jetTrkMax/hadd.log"
