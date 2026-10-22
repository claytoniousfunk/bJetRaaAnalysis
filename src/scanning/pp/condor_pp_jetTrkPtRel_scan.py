import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# pp track ptRel for the data-driven light template (all tracks and the leading
# one), 2026-10-19. The macro picks its own input list (the HighEGJet forests of
# pp_scan.C) and output directory; the submit script asks it with
# pp_jetTrkPtRel_scan.C(0). Run from src/scanning/pp/.
# Merge afterwards:  hadd -f <outdir>.root <outdir>/pp_scan_output_*.root
exe = 'pp_jetTrkPtRel_scan.C'
jobname = None   # default: <macro>_<list name>

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 10,               # forest files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = False        # trial-run script_0.sh first, then condor_submit
)
