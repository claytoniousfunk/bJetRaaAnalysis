import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# Same scheme as condor_PYTHIAHYDJET_scan.py: the input list is whatever
# PYTHIAHYDJET_scan_response.C reads (asked via its query mode), so the job count
# always matches the list the jobs will open. The output directory follows the
# flags in config_PYTHIAHYDJET.h plus CENT_SCHEME_SUFFIX from config_centrality.h
# (e.g. _ultraFineCentBins) and must exist on EOS before submitting.
#
# Run from src/scanning/PYTHIAHYDJET/.
exe = 'PYTHIAHYDJET_scan_response.C'
jobname = None   # default: <macro>_<list name>; set a string to override

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 5,                # input files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = False        # trial-run script_0.sh first; this scan books ~4x more histograms under ultraFine
)
