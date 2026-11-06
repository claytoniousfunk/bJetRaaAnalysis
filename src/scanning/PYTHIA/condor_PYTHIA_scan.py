import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# input list is whatever the macro's sample flags select (set them in its config
# header); the submit script asks the macro, so there is nothing to match here
exe = 'PYTHIA_scan.C'
jobname = None   # default: <macro>_<list name>; set a string to override

# Event parity (2026-11-06): `python3 condor_PYTHIA_scan.py even` (or `odd`) scans only
# the even- (odd-) indexed events, PYTHIA_scan's parity argument, for independent
# template / pseudo-data halves. Output names and job directories carry _evenEvents /
# _oddEvents. No argument = every event, as before.
parity = sys.argv[1] if len(sys.argv) > 1 else ''
if parity not in ('', 'even', 'odd'):
    raise SystemExit("usage: condor_PYTHIA_scan.py [even|odd]")
macro_args = {'': '', 'even': ', 0', 'odd': ', 1'}[parity]
if parity and jobname is None:
    jobname = f'PYTHIA_scan_{parity}Events'

submitScan(exe, jobname,
	time_flavour = 'workday',  # 8h
	nsplit = 5,                # input files per condor job
	njobs_max = None,          # e.g. 1 for a trial run
	auto_submit = True,
	macro_args = macro_args
)
