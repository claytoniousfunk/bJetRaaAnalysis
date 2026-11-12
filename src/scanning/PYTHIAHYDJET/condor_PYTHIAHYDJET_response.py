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
#
# Event parity (2026-10-08): `python3 condor_PYTHIAHYDJET_response.py even` (or
# `odd`) runs only that half -- the macro's parity argument, overriding
# onlyEvenEvents / onlyOddEvents in config_PYTHIAHYDJET.h -- with its own job
# directory, so both halves can be queued together. No argument: the config's
# setting (normally all events).
#
# Sample (2026-11-11): name one or more of DiJet / MuJet / BJet to queue those
# samples in one go, each with its own job directory and output name
# (PYTHIAHYDJET_response_<sample>_...), e.g.
#     python3 condor_PYTHIAHYDJET_response.py MuJet BJet
#     python3 condor_PYTHIAHYDJET_response.py even MuJet BJet
#     python3 condor_PYTHIAHYDJET_response.py even odd DiJet MuJet BJet
# Naming both even and odd queues both halves of every sample (one job
# directory per sample x half); the halves split on the event index within each
# forest file, so even + odd is the full sample.
# The sample is passed to every job (the macro's 3rd argument), so the config's
# sample flags are ignored and can change while the jobs are queued; every other
# config flag is read when each job starts, so leave config_PYTHIAHYDJET.h alone
# until they finish. No sample named: the config's sample flag, as before.
exe = 'PYTHIAHYDJET_scan_response.C'

parities = {'even': 0, 'odd': 1}
sampleIdx = {'DiJet': 1, 'MuJet': 2, 'BJet': 3}
parityList, samples = [], []
for a in sys.argv[1:]:
    if a in parities and a not in parityList: parityList.append(a)
    elif a in sampleIdx and a not in samples: samples.append(a)
    else: raise SystemExit("usage: condor_PYTHIAHYDJET_response.py [even] [odd] [DiJet] [MuJet] [BJet]")
if not parityList: parityList = ['']  # config's setting

def submit(macro_args, jobname):
    submitScan(exe, jobname,
        time_flavour = 'workday',  # 8h
        nsplit = 5,                # input files per condor job
        njobs_max = None,          # e.g. 1 for a trial run
        auto_submit = True,        # submits on run; set False to trial-run script_0.sh first (this scan books ~4x more histograms under ultraFine)
        macro_args = macro_args
    )

for parity in parityList:
    parTag = '_' + parity + 'Events' if parity else ''
    if not samples:
        # jobname None = <macro>_<list name>
        submit({'': '', 'even': ', 0', 'odd': ', 1'}[parity],
               os.path.splitext(exe)[0] + parTag if parity else None)
    else:
        for name in samples:
            submit(f', {parities.get(parity, -1)}, {sampleIdx[name]}',
                   f'{os.path.splitext(exe)[0]}_{name}{parTag}')
