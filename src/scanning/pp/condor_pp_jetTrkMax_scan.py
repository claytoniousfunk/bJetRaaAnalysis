#!/usr/bin/env python3
"""Submit pp_jetTrkMax_scan.C to condor, one job per block of skim groups.

Replaces run_pp_jetTrkMax_scan.sh (an interactive xargs loop on lxplus, which
dies with the ssh session: the 2026-10-09 run left 5 group outputs). 2026-10-19.

The macro takes a skim group number and builds its input and output paths
itself (pp_jetTrkMax_scan.C, lines ~96-97), so this script only needs the list
of groups, read from the skim directory at submit time. Each job runs `ngroup`
groups in turn, interpreted as the 2025 scan was (ACLiC would lack the implicit
`using namespace std` that eventMap.h relies on), so there is no shared .so for
concurrent jobs to corrupt.

The macro does not check its input, so a failed group writes no file and ROOT
may still exit 0. Each job deletes a group's old output before running it and
exits 1 if it is not there afterwards, so condor reports the job as failed.

Run from src/scanning/pp/ (the job scripts cd here; the include paths are
relative to it). Merge afterwards with the hadd command printed at the end.
"""

import math
import os
import re

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

jobname = 'pp_jetTrkMax_trkPtAxis'

skimdir = '/eos/cms/store/group/phys_heavyions/cbennett/skims/output_skims_pp_HighEGJet_withJetTrackMaxInfo'

ngroup = 20                   # skim groups per job; ~77k events, ~3 MB each
time_flavour = '"longlunch"'  # 2h
request_memory = 2000         # MB

# Same view as the other pp/PbPb condor scripts and ~/.bashrc.
env_setup = 'source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh'

auto_submit = False           # run the trial job printed at the end first
njobs_max = None              # e.g. 1 for a trial; None = all groups

# ---------------------------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))
macro = 'pp_jetTrkMax_scan.C'
macro_path = os.path.join(here, macro)
if not os.path.isfile(macro_path):
    raise SystemExit(f'ERROR: {macro} not found in {here}')

# the macro's own output directory, so the two can never disagree
with open(macro_path) as f:
    m = re.search(r'"([^"]*output_pp_jetTrkMax[^"/]*)/pp_scan_output_', f.read())
if not m:
    raise SystemExit(f'ERROR: no output directory found in {macro}')
outdir = m.group(1)

if not os.path.isdir(skimdir):
    raise SystemExit(f'ERROR: skim directory not found: {skimdir}')
groups = sorted(int(g.group(1)) for g in
                (re.match(r'^pp_skim_output_(\d+)\.root$', n) for n in os.listdir(skimdir)) if g)
if not groups:
    raise SystemExit(f'ERROR: no pp_skim_output_<n>.root in {skimdir}')

njobs = int(math.ceil(float(len(groups)) / ngroup))
if njobs_max is not None:
    njobs = min(njobs, njobs_max)

jobdir = os.path.join(here, jobname)
os.makedirs(os.path.join(jobdir, 'log'), exist_ok=True)
os.makedirs(outdir, exist_ok=True)

for i in range(njobs):
    block = ' '.join(str(g) for g in groups[i * ngroup:(i + 1) * ngroup])
    lines = ['#!/bin/bash', '', env_setup, '', f'cd {here}', '',
             'fail=0',
             f'for g in {block}; do',
             f'  out={outdir}/pp_scan_output_$g.root',
             '  rm -f "$out"',
             f"  root -l -b -q \"{macro}($g)\"",
             '  if [ ! -s "$out" ]; then echo "ERROR: group $g wrote no output"; fail=1; fi',
             'done',
             'exit $fail',
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
+JobFlavour = {time_flavour}
Queue {njobs}
"""
sub_path = os.path.join(jobdir, 'condor_submit.cfg')
with open(sub_path, 'w') as f:
    f.write(sub)

merged = f'{os.path.dirname(outdir)}/{os.path.basename(outdir)}.root'
print(f'Generated {njobs} jobs for {len(groups)} skim groups (ngroup={ngroup}), '
      f'groups {groups[0]}-{groups[-1]}')
print(f'Output: {outdir}/pp_scan_output_<group>.root')
print(f'Submit file: {sub_path}')
if auto_submit:
    os.system(f'condor_submit {sub_path}')
else:
    print('\nauto_submit is False. Trial-run one job locally first (time it):')
    print(f'  time {jobdir}/script_0.sh')
    print('then submit with:')
    print(f'  condor_submit {sub_path}')
print('\nWhen all jobs are done, check the count and merge:')
print(f'  ls {outdir}/pp_scan_output_*.root | wc -l    # expect {len(groups)}')
print(f'  hadd -f {merged} {outdir}/pp_scan_output_*.root')
