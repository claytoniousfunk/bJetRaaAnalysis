#!/usr/bin/env python3
"""Submit caloVsPFJetScale_scan.C to condor, for PbPb or pp.

Adapted from src/scanning/PbPb/condor_PbPb_muonTagAndProbe.py; the job
mechanics are the same. Set `system` below and run once per system:

    system = 'PbPb'   HIHardProbes, akPu4Calo vs akCs4PF
    system = 'pp'     HighEGJet,    ak4Calo   vs ak4PF

pp uses the list pp_scan.C reads for HighEGJet. PbPb uses the _fresh
HardProbes forest, a re-forest of the dataset PbPb_scan.C was scanned on
(its list still points at the deleted original).

Run from src/scanning/jetCollection/ (the job scripts cd here, and the
macro's JEC paths are relative to it).

Merge afterwards into the names the plot macro reads:
    hadd -f PbPb_caloVsPFJetScale.root            <PbPb outdir>/cvpf_*.root
    hadd -f pp_HighEGJet_caloVsPFJetScale.root    <pp outdir>/cvpf_*.root
copy them to rootFiles/scanningOuput/PbPb/ and .../pp/, and run
src/plots/jetPt/jetCollection/plotCaloOverPFJetScale_PbPbVsPP.C.
provenance is not mergeable; hadd keeps the first copy, which is fine since
every job writes the same one apart from the input list.

Statistics: the quantity that matters is 50-80% above 200 GeV. The full
PbPb calo Jet100 scan (2026-09-20) has ~11k calo jets at 200-320 GeV in
50-80% (~0.9k at 320-500), so with a ratio RMS of ~0.1 the mean is known to
~0.1%. A tenth of the files (njobs_max) is already enough for a first look.
"""

import math
import os

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

system = 'PbPb'               # 'PbPb' or 'pp'

if system == 'PbPb':
    # the _fresh forest (crab 260926_022130); the original
    # fileNames_HIHardProbes_withCaloAndFlowJets.txt points at files since deleted
    dblist = '../../../fileNames/fileNames_HIHardProbes_withCaloAndFlowJets_fresh_partial.txt'
    dataset = 'HardProbes'
    isPP = False
elif system == 'pp':
    dblist = '../../../fileNames/fileNames_pp_HighEGJet.txt'
    dataset = 'HighEGJet'
    isPP = True
else:
    raise SystemExit(f"ERROR: system must be 'PbPb' or 'pp', not {system!r}")

# <system>_<dataset>_<study>, as in the muonTagAndProbe scripts
jobname = f'{system}_{dataset}_caloVsPFJetScale'

# Must be a directory the jobs can write: EOS via the FUSE mount works from
# lxplus condor.
outdir = f'/eos/cms/store/group/phys_heavyions/cbennett/scanningOutput/output_{jobname}'

nsplit = 50                   # input files per job; time the trial and adjust
time_flavour = '"workday"'    # 8h
request_memory = 2000         # MB

# Same view as the other pp/PbPb condor scripts and ~/.bashrc. A mismatch puts
# two ROOT installations in one process.
env_setup = 'source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh'

auto_submit = False           # run the trial job printed at the end first
njobs_max = None              # e.g. 1 for a trial; None = all files

# ---------------------------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))
wrapper = 'run_caloVsPFJetScale_condor.C'
for need in (wrapper, 'caloVsPFJetScale_scan.C'):
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

    out = os.path.join(outdir, f'cvpf_{i}.root')
    lines = ['#!/bin/bash', 'set -e', '', env_setup, '', f'cd {here}', '',
             f"root -l -b -q '{wrapper}(\"{sub_list}\", \"{out}\", {'true' if isPP else 'false'})'",
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

print(f'{system}: generated {njobs} jobs for {nfiles} input files (nsplit={nsplit})')
print(f'Output: {outdir}/cvpf_<i>.root')
print(f'Submit file: {sub_path}')
if auto_submit:
    os.system(f'condor_submit {sub_path}')
else:
    print('\nauto_submit is False. Trial-run one job locally first (time it):')
    print(f'  time {jobdir}/script_0.sh')
    print('then submit with:')
    print(f'  condor_submit {sub_path}')
