import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# PYTHIA+HYDJET template scans with pThat > 15, one per sample, all queued at
# once; the pThat > 15 counterpart of condor_PYTHIAHYDJET_scan_pThat30_series.py.
# Sample and pThat cut are passed to the macro per job (PYTHIAHYDJET_scan
# arguments 2 and 3), so config_PYTHIAHYDJET.h's sample flags are ignored and
# can be changed while these jobs are in the queue. Every other setting is still
# read from the config when each job starts: do not edit it until they finish.
# The output names carry _pThat-15, so they do not overwrite the pThat-30 scans.
# The DiJet job writes to the same directory as condor_PYTHIAHYDJET_scan.py with
# the DiJet config flag (both are pThat-15 DiJet): do not run both on one day.
exe = 'PYTHIAHYDJET_scan.C'
pThatMin = 15.
samples = {'DiJet': 1, 'MuJet': 2, 'BJet': 3}

for name, idx in samples.items():
    submitScan(exe, f'PYTHIAHYDJET_scan_{name}_pThat-{pThatMin:.0f}',
        time_flavour = 'workday',  # 8h
        nsplit = 5,                # input files per condor job
        njobs_max = None,          # e.g. 1 for a trial run
        auto_submit = True,
        macro_args = f', {idx}, {pThatMin}'
    )
