#!/usr/bin/env bash
# Run pp_jetTrkMax_scan.C over every legacy skim group, in parallel, then hadd.
# Run ON LXPLUS from src/scanning/pp/ (2026-10-09). The skims are small (~77k
# events, ~3 MB each), so this is a local loop rather than condor.
#   ./run_pp_jetTrkMax_scan.sh [nParallel=8] [firstN=all]   # firstN: trial on the first N groups
set -euo pipefail
NP=${1:-8}; FIRST=${2:-0}
SKIMS=/eos/cms/store/group/phys_heavyions/cbennett/skims/output_skims_pp_HighEGJet_withJetTrackMaxInfo
OUT=$(grep -o '/eos/[^"]*output_pp_jetTrkMax[^"/]*' pp_jetTrkMax_scan.C | head -1)   # the macro's own output dir
mkdir -p "$OUT" logs_jetTrkMax
GROUPS=$(ls $SKIMS | sed -n 's/^pp_skim_output_\([0-9]*\)\.root$/\1/p' | sort -n)
if [ "$FIRST" -gt 0 ]; then GROUPS=$(echo "$GROUPS" | head -n "$FIRST"); fi
echo "groups: $(echo "$GROUPS" | wc -l), parallel $NP, output $OUT"
# interpreted, as the 2025 scan (ACLiC would lack the implicit `using namespace std`
# that eventMap.h relies on)
echo "$GROUPS" | xargs -P "$NP" -I{} sh -c "root -l -b -q 'pp_jetTrkMax_scan.C({})' > logs_jetTrkMax/g{}.log 2>&1 || echo 'group {} FAILED'"
echo "done: $(ls $OUT/pp_scan_output_*.root 2>/dev/null | wc -l) outputs; failures: $(grep -L 'evt frac\|Processing events' logs_jetTrkMax/g*.log | wc -l)"
hadd -f "$OUT/../$(basename $OUT).root" $OUT/pp_scan_output_*.root > logs_jetTrkMax/hadd.log && echo "merged: $OUT/../$(basename $OUT).root"
