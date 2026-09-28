#!/usr/bin/env python3
"""Submit PbPb_caloTowerAnalyzer.C to condor, one input file per job.

The CaloTower counterpart of condor_PbPb_pfCandAnalyzer.py, and deliberately a
near-copy of it: same standalone-binary default, same env_setup, same Delphes
note, same trial-run workflow. Differences are only where the two analyzers
actually differ:

  * it builds and runs caloTowerAnalyzer_PbPb (Makefile.caloTowerAnalyzer),
    which is a separate target from pfCandAnalyzer_PbPb -- the two share every
    analysis header, so one target name would let each silently overwrite the
    other's binary;
  * the job scripts do NOT pass a file-list override, so the analyzer picks its
    list from config_PbPb.h exactly as in production. The override exists for
    local smoke tests (see the module docstring in PbPb_caloTowerAnalyzer.C);
  * request_memory is lower, see the note below.

Run from src/scanning/PbPb/ (the analyzer resolves its input list relative to
that directory).

BEFORE SUBMITTING, check that headers/AnalysisSetup/caloTowers.h says what you
mean. towerEtMin and towerEtaMaxCluster change which towers are clustered -- 76%
of towers clear the default 0.3 GeV and the |eta| < 3.0 cut removes HF -- and
they go into the output filename, so two settings cannot be merged by accident.
"""

import math
import os

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

jobname = 'PbPb_caloTowerAnalyzer'

# Input list -- MUST be the one PbPb_caloTowerAnalyzer.C selects from
# config_PbPb.h, because the jobs pass a line index that is matched against
# `ifile` inside the analyzer. A mismatch here silently scans different files
# than the indices claim.
#
# This is the only forest production that carries rechitanalyzerpp towers:
# forest_HIMinimumBias0_Part1_withCaloAndFlowJets_withPFAndTowers_fresh, CRAB
# 260926, 1992 files. It requires doMinBiasSample_Part1 = true in config_PbPb.h
# with every other sample flag false; the analyzer refuses to run otherwise.
dblist = '../../../fileNames/fileNames_HIMinimumBias0_Part1_withCaloAndFlowJets_withPFAndTowers_fresh.txt'

# How the analyzer is run.
#
#   True  (default): build a standalone binary once here, jobs run it directly.
#                    No cling, so ROOT never autoloads Delphes and its embedded
#                    FastJet cannot hijack fastjet::ClusterSequence (that
#                    collision segfaults the ACLiC path). Also removes ~2000
#                    redundant ACLiC compiles.
#   False: legacy path -- each job compiles the macro with ACLiC via
#          run_caloTowerAnalyzer_condor.C. KNOWN BROKEN on any LCG view that
#          ships Delphes: it segfaults in ClusterSequence's constructor. Kept
#          only for use on a stack without Delphes.
use_standalone_binary = True

binary = 'caloTowerAnalyzer_PbPb'
makefile = 'Makefile.caloTowerAnalyzer'
exe = 'run_caloTowerAnalyzer_condor.C'   # used only when use_standalone_binary is False

# Sourced at the top of every job script. Must put fastjet-config on PATH and
# libfastjet on LD_LIBRARY_PATH. Replace with the view you use interactively --
# `which fastjet-config` in a working shell will tell you which one that is.
# Must match the view the binary was built against -- LCG_106 (ROOT 6.32.02),
# which is what ~/.bashrc sources. A mismatched view puts two ROOT installations
# in one process and segfaults in TCling's constructor. See the HYDJET copy of
# this script for the full note.
env_setup = 'source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh'

# The LCG views ship Delphes, which bundles its own FastJet build and claims the
# fastjet:: namespace in its rootmap. ROOT's class autoloading can therefore
# resolve fastjet::ClusterSequence into libDelphesDisplay.so instead of the
# standalone libfastjet the macro was compiled against, which segfaults.
# run_caloTowerAnalyzer_condor.C preloads the real libfastjet to prevent this.
# Set this True to also force it at the process level, which is airtight.
force_fastjet_preload = False

time_flavour = '"tomorrow"'   # 1 day; mixed-event running is slow, workday (8h) may not be enough

# MB. The mixed-event pool holds N_mixedEventsInPool events' worth of CLUSTERED
# towers -- the pool is filtered at fill time, so it carries ~1200 of the ~2250
# towers in an event rather than all of them. That is smaller than the PF
# candidate pool, but the measurement behind it came from a forest that is
# mostly peripheral, and central events have considerably more towers. 3000 is
# kept as the pfCandAnalyzer value rather than tuned down on one file's worth of
# evidence; lower it once a trial job's real usage is known.
request_memory = 3000

nsplit = 1                    # input files per condor job

# Full submission: all 1992 files. Left at auto_submit = False so the trial job
# printed at the end is run first -- this is a fresh forest production and the
# first job is the only cheap chance to find out it is not what we think.
auto_submit = False           # set True once a trial job has succeeded
njobs_max = None              # e.g. 1 for a trial run; None = all files

# ---------------------------------------------------------------------------

here = os.path.dirname(os.path.abspath(__file__))

dblist_path = dblist if os.path.isabs(dblist) else os.path.join(here, dblist)
if not os.path.isfile(dblist_path):
    raise SystemExit(f'ERROR: input list not found: {dblist_path}')
if use_standalone_binary:
    # Build once, now, in the submitting shell -- so a compile failure is seen
    # here instead of 1993 times on the farm.
    print(f'Building {binary} ...')
    rc = os.system(f'cd {here} && make -f {makefile}')
    if rc != 0:
        raise SystemExit(f'ERROR: build failed (make -f {makefile} returned {rc})')
    if not os.path.isfile(os.path.join(here, binary)):
        raise SystemExit(f'ERROR: build reported success but {binary} is missing')
    print(f'Built {os.path.join(here, binary)}\n')
elif not os.path.isfile(os.path.join(here, exe)):
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

    lines = ['#!/bin/bash', 'set -e', '', env_setup, '']
    if force_fastjet_preload and not use_standalone_binary:
        lines += ['export LD_PRELOAD="$(fastjet-config --prefix)/lib/libfastjet.so"', '']
    lines += [f'cd {here}', '']
    for idx in range(start, end + 1):
        if use_standalone_binary:
            lines.append(f'./{binary} {idx}')
        else:
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
# Deliberately False: env_setup is sourced at the top of every job script and is
# meant to be the whole environment. With getenv = True the job inherits the
# submitting shell too and layers the view on top of it.
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
