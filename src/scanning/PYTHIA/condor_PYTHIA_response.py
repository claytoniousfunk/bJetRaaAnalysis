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
#
# Event parity (2026-11-11, as condor_PYTHIAHYDJET_response.py):
#     python3 condor_PYTHIA_response.py            all events (config's onlyEvenEvents / onlyOddEvents)
#     python3 condor_PYTHIA_response.py even       even events only
#     python3 condor_PYTHIA_response.py even odd   both halves, each with its own job directory
# The parity is the macro's 2nd argument, overriding onlyEvenEvents / onlyOddEvents in
# config_PYTHIA.h; the output names carry _evenEvents / _oddEvents. The halves split on the
# event index within each forest file, so even + odd is the full sample. Every other config flag
# is read when each job starts, so leave config_PYTHIA.h alone until the jobs finish.
exe = 'PYTHIA_scan_response.C'

parities = {'even': 0, 'odd': 1}
parityList = []
for a in sys.argv[1:]:
    if a in parities and a not in parityList: parityList.append(a)
    else: raise SystemExit("usage: condor_PYTHIA_response.py [even] [odd]")
if not parityList: parityList = ['']  # config's setting

for parity in parityList:
    submitScan(exe,
        os.path.splitext(exe)[0] + '_' + parity + 'Events' if parity else None,   # None = <macro>_<list name>
        time_flavour = 'workday',  # 8h
        nsplit = 5,                # input files per condor job
        njobs_max = None,          # e.g. 1 for a trial run
        auto_submit = True,
        macro_args = {'': '', 'even': ', 0', 'odd': ', 1'}[parity]
    )
