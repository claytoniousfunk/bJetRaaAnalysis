import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# PbPb leading-track ptRel (data-driven light template), 2026-10-19. The macro picks
# its own input list (the PbPb_pfCandAnalyzer.C MinBias forests with tracks; the
# leading track is found from PbPbTracks/trackTree) and output directory; the
# submit script asks it with PbPb_jetTrkMax_scan.C(0). Run from src/scanning/PbPb/.
# Merge afterwards:  hadd -f <outdir>.root <outdir>/PbPb_scan_output_*.root
exe = 'PbPb_jetTrkMax_scan.C'
jobname = None   # default: <macro>_<list name>

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 10,               # forest files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = False        # trial-run script_0.sh first, then condor_submit
)
