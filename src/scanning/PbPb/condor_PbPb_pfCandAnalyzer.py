#!/usr/bin/env python3
"""Submit PbPb_pfCandAnalyzer.C to condor, one input file per job.

Modelled on src/scanning/HYDJET/condor_HYDJET_pfCandAnalyzer.py, with three
changes for the PbPb mixed-event/FastJet running:

  * the environment setup is written into each job script instead of relying on
    `getenv = True`, so a job cannot silently inherit a shell without FastJet;
  * it drives run_pfCandAnalyzer_condor.C, which uses a private ACLiC build
    directory (2000 jobs building in one directory corrupt each other's .so);
  * auto_submit defaults to False, and njobs_max lets you queue a handful of
    jobs first. Do a single-job trial before queueing the full list.

Run from src/scanning/PbPb/ (the analyzer resolves its input list relative to
that directory).
"""

import math
import os

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

jobname = 'PbPb_pfCandAnalyzer_mixedEvent_fastJet'

# input list -- must be the same one PbPb_pfCandAnalyzer.C selects via its
# doMinBiasSample / doHardProbesSample / doSingleMuonSample flags, because the
# job index is matched against `ifile` inside the analyzer.
dblist = '../../../fileNames/fileNames_HIMinimumBias0_Part1_withTracksAndPFCandidates.txt'

exe = 'run_pfCandAnalyzer_condor.C'

# Sourced at the top of every job script. Must put fastjet-config on PATH and
# libfastjet on LD_LIBRARY_PATH. Replace with the view you use interactively --
# `which fastjet-config` in a working shell will tell you which one that is.
env_setup = 'source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh'

time_flavour = '"tomorrow"'   # 1 day; mixed-event running is slow, workday (8h) may not be enough
request_memory = 3000         # MB; the 100-event PF candidate pool is ~30 MB, ROOT buffers dominate
nsplit = 1                    # input files per condor job

auto_submit = False           # set True once a trial job has succeeded
njobs_max = None              # e.g. 1 for a trial run; None = all files

# ---------------------------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))

dblist_path = dblist if os.path.isabs(dblist) else os.path.join(here, dblist)
if not os.path.isfile(dblist_path):
    raise SystemExit(f'ERROR: input list not found: {dblist_path}')
if not os.path.isfile(os.path.join(here, exe)):
    raise SystemExit(f'ERROR: executable macro not found: {os.path.join(here, exe)}')

with open(dblist_path) as f:
    files = [l.strip() for l in f if l.strip()]

nfiles = len(files)
njobs = int(math.ceil(float(nfiles) / nsplit))
if njobs_max is not None:
    njobs = min(njobs, njobs_max)

jobdir = os.path.join(here, jobname)
os.makedirs(os.path.join(jobdir, 'log'), exist_ok=True)

# one script per job
for i in range(njobs):
    start = i * nsplit + 1
    end = min((i + 1) * nsplit, nfiles)

    lines = ['#!/bin/bash', 'set -e', '', env_setup, '', f'cd {here}', '']
    for idx in range(start, end + 1):
        lines.append(f"root -l -b -q '{exe}({idx})'")
    lines.append('')

    script_path = os.path.join(jobdir, f'script_{i}.sh')
    with open(script_path, 'w') as f:
        f.write('\n'.join(lines))
    os.chmod(script_path, 0o755)

sub = f"""Universe   = vanilla
Executable = {jobdir}/script_$(ProcId).sh
Log        = {jobdir}/log/job_$(ProcId).log
Output     = {jobdir}/log/job_$(ProcId).out
Error      = {jobdir}/log/job_$(ProcId).err
getenv     = True
request_memory = {request_memory}
x509userproxy = $ENV(X509_USER_PROXY)
use_x509userproxy = True
+JobFlavour = {time_flavour}
Queue {njobs}
"""

sub_path = os.path.join(jobdir, 'condor_submit.cfg')
with open(sub_path, 'w') as f:
    f.write(sub)

print(f'Generated {njobs} jobs for {nfiles} input files (nsplit={nsplit})')
if njobs_max is not None:
    print(f'NOTE: njobs_max={njobs_max} -- only the first {njobs} of '
          f'{int(math.ceil(float(nfiles) / nsplit))} jobs were generated')
print(f'Submit file: {sub_path}')

if auto_submit:
    os.system(f'condor_submit {sub_path}')
else:
    print('\nauto_submit is False. Trial-run one job locally first:')
    print(f'  {jobdir}/script_0.sh')
    print('then submit with:')
    print(f'  condor_submit {sub_path}')
