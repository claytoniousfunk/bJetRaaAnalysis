import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# Same scheme as condor_PYTHIA_scan.py: the input list and output directory are
# whatever PYTHIA_scan_response.C reports in query mode (its config_PYTHIA.h
# flags), so the job count always matches the list the jobs open and the output
# directory is created before submitting.
#
# Run from src/scanning/PYTHIA/.
exe = 'PYTHIA_scan_response.C'
jobname = None   # default: <macro>_<list name>; set a string to override

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 5,                # input files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = True
)
