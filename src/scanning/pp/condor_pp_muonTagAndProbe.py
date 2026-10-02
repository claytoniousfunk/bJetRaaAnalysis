#!/usr/bin/env python3
"""Submit pp_muonTagAndProbe_scan.C to condor.

Unlike the production scans (scan_condor.py), the tag-and-probe scan takes its
file list and output path as arguments rather than choosing them from a config
header. So this script splits the list itself: each job gets its own sub-list
of `nsplit` files under <jobname>/lists/ and writes one output file to outdir.

Run from src/scanning/pp/ (the job scripts cd here, and the macro's include
paths are relative to it).

Each job compiles the macro with ACLiC into its own scratch directory via
run_muonTagAndProbe_condor.C (~30 s), which keeps concurrent jobs from
corrupting a shared .so.

Merge afterwards with
    hadd -f pp_HighEGJet_tnp.root <outdir>/tnp_*.root
The histograms and the tnp tree add; provenance is not mergeable and hadd keeps
the first copy, which is fine since every job writes the same one apart from the
input list.
"""

import math
import os

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

jobname = 'pp_HighEGJet_muonTagAndProbe'

# HighEGJet forest: photon/jet triggers, so the probe muon is unbiased by the
# trigger. Do NOT point this at SingleMuon -- see the scan's header.
dblist = '../../../fileNames/fileNames_pp_HighEGJet.txt'
isMC = False

# Must be a directory the jobs can write: EOS via the FUSE mount works from
# lxplus condor.
outdir = '/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/output_pp_HighEGJet_muonTagAndProbe'

nsplit = 50                   # input files per job; ~70k events per file
time_flavour = '"workday"'    # 8h
request_memory = 2000         # MB; the scan holds no event pool

# Same view as the other pp/PbPb condor scripts and ~/.bashrc. A mismatch puts
# two ROOT installations in one process.
env_setup = 'source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh'

auto_submit = False           # run the trial job printed at the end first
njobs_max = None              # e.g. 1 for a trial; None = all files

# ---------------------------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))
wrapper = 'run_muonTagAndProbe_condor.C'
for need in (wrapper, 'pp_muonTagAndProbe_scan.C'):
    if not os.path.isfile(os.path.join(here, need)):
        raise SystemExit(f'ERROR: {need} not found in {here}')

dblist_path = dblist if os.path.isabs(dblist) else os.path.join(here, dblist)
if not os.path.isfile(dblist_path):
    raise SystemExit(f'ERROR: input list not found: {dblist_path}')
with open(dblist_path) as f:
    files = [l.strip() for l in f if l.strip()]

nfiles = len(files)
njobs = int(math.ceil(float(nfiles) / nsplit))
if njobs_max is not None:
    njobs = min(njobs, njobs_max)

jobdir = os.path.join(here, jobname)
os.makedirs(os.path.join(jobdir, 'log'), exist_ok=True)
os.makedirs(os.path.join(jobdir, 'lists'), exist_ok=True)
os.makedirs(outdir, exist_ok=True)

for i in range(njobs):
    sub_list = os.path.join(jobdir, 'lists', f'list_{i}.txt')
    with open(sub_list, 'w') as f:
        f.write('\n'.join(files[i * nsplit:(i + 1) * nsplit]) + '\n')

    out = os.path.join(outdir, f'tnp_{i}.root')
    lines = ['#!/bin/bash', 'set -e', '', env_setup, '', f'cd {here}', '',
             f"root -l -b -q '{wrapper}(\"{sub_list}\", \"{out}\", {'true' if isMC else 'false'})'",
             '']
    script_path = os.path.join(jobdir, f'script_{i}.sh')
    with open(script_path, 'w') as f:
        f.write('\n'.join(lines))
    os.chmod(script_path, 0o755)

sub = f"""Universe   = vanilla
Executable = {jobdir}/script_$(ProcId).sh
Log        = {jobdir}/log/job_$(ProcId).log
Output     = {jobdir}/log/job_$(ProcId).out
Error      = {jobdir}/log/job_$(ProcId).err
# env_setup in each script is the whole environment; do not inherit the shell's
getenv     = False
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
print(f'Output: {outdir}/tnp_<i>.root')
print(f'Submit file: {sub_path}')
if auto_submit:
    os.system(f'condor_submit {sub_path}')
else:
    print('\nauto_submit is False. Trial-run one job locally first (time it):')
    print(f'  time {jobdir}/script_0.sh')
    print('then submit with:')
    print(f'  condor_submit {sub_path}')
