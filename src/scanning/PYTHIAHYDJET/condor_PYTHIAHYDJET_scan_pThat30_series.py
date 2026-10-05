import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..')))
from myProcesses.condor.scan_condor import submitScan

# PYTHIA+HYDJET template scans with pThat > 30, one per sample, all queued at
# once. Sample and pThat cut are passed to the macro per job (PYTHIAHYDJET_scan
# arguments 2 and 3), so config_PYTHIAHYDJET.h's sample flags are ignored and
# can be changed while these jobs are in the queue. Every other setting is still
# read from the config when each job starts: do not edit it until they finish.
# The output names carry _pThat-30, so they do not overwrite the pThat-15 scans.
exe = 'PYTHIAHYDJET_scan.C'
pThatMin = 30.
samples = {'DiJet': 1, 'MuJet': 2, 'BJet': 3}

for name, idx in samples.items():
    submitScan(exe, f'PYTHIAHYDJET_scan_{name}_pThat-{pThatMin:.0f}',
        time_flavour = 'workday',  # 8h
        nsplit = 5,                # input files per condor job
        njobs_max = None,          # e.g. 1 for a trial run
        auto_submit = True,
        macro_args = f', {idx}, {pThatMin}'
    )
