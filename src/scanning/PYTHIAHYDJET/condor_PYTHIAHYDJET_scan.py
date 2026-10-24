import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# input list is whatever the macro's sample flags select (set them in its config
# header); the submit script asks the macro, so there is nothing to match here
exe = 'PYTHIAHYDJET_scan.C'
jobname = None   # default: <macro>_<list name>; set a string to override
# pThat cut passed to each job (PYTHIAHYDJET_scan argument 3); None = pthatcut from
# common.h (15). 30 matches the template scans (2026-10-24, CSVv2 DiJet scan). The
# output name carries _pThat-30, so it does not overwrite the pThat-15 scans.
pThatMin = 30.

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 5,                # input files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = True,
	macro_args = '' if pThatMin is None else f', 0, {pThatMin}'   # sample 0 = config flags
)
