"""Condor submission for the plain (non-FastJet) scan macros.

Used by condor_{PYTHIA,pp,PbPb,PYTHIAHYDJET}_scan.py. Each job runs
`root -l -b -q 'macro.C(idx)'` for a block of `nsplit` consecutive 1-based file
indices. The macros pick their own input list (from the sample flags in their
config headers) and write to EOS themselves. The list is obtained by running the
macro in query mode (`macro.C(0)`), so the job count always matches what the
jobs will read; set the sample in the macro's config header, nowhere else.

Run from the macro's directory: the macros resolve their list relative to it.
"""

import math
import os
import re
import subprocess


def _queryList(exe):
    """Ask the macro which input list its current flags select (exe(0) prints it)."""
    out = subprocess.run(['root', '-l', '-b', '-q', f'{exe}(0)'],
                         capture_output=True, text=True).stdout
    m = re.search(r'INPUTFILELIST=(\S+)', out)
    if not m or not m.group(1):
        raise SystemExit(f'ERROR: {exe}(0) did not report an input list. Output tail:\n'
                         + '\n'.join(out.splitlines()[-15:]))
    return m.group(1)


def submitScan(exe, jobname=None, time_flavour='workday', nsplit=5,
               njobs_max=None, auto_submit=True):
    pwd = os.getcwd()
    if not os.path.isfile(exe):
        raise SystemExit(f'ERROR: macro not found: {os.path.join(pwd, exe)}')

    dblist = _queryList(exe)
    if not os.path.isfile(dblist):
        raise SystemExit(f'ERROR: input list not found: {os.path.abspath(dblist)}')
    print(f'{exe} reads {dblist}')
    if jobname is None:
        jobname = os.path.splitext(exe)[0] + '_' + os.path.splitext(os.path.basename(dblist))[0]

    with open(dblist) as f:
        nfiles = sum(1 for l in f if l.strip())
    njobs = int(math.ceil(float(nfiles) / nsplit))
    if njobs_max is not None:
        njobs = min(njobs, njobs_max)

    jobdir = f'{pwd}/{jobname}'
    os.makedirs(f'{jobdir}/log', exist_ok=True)

    for i in range(njobs):
        start = i * nsplit + 1
        end = min((i + 1) * nsplit, nfiles)
        script = f'#!/bin/bash\ncd {pwd}\n'
        for idx in range(start, end + 1):
            script += f"root -l -b -q '{exe}({idx})'\n"
        path = f'{jobdir}/script_{i}.sh'
        with open(path, 'w') as f:
            f.write(script)
        os.chmod(path, 0o755)

    sub = f"""Universe   = vanilla
Executable = {jobdir}/script_$(ProcId).sh
Log        = {jobdir}/log/job_$(ProcId).log
Output     = {jobdir}/log/job_$(ProcId).out
Error      = {jobdir}/log/job_$(ProcId).err
getenv     = True
x509userproxy = $ENV(X509_USER_PROXY)
use_x509userproxy = True
+JobFlavour = "{time_flavour}"
Queue {njobs}
"""
    sub_path = f'{jobdir}/condor_submit.cfg'
    with open(sub_path, 'w') as f:
        f.write(sub)

    print(f'Generated {njobs} jobs for {nfiles} input files (nsplit={nsplit})')
    if njobs_max is not None:
        print(f'NOTE: njobs_max={njobs_max} -- jobs capped')
    print(f'Submit file: {sub_path}')
    if auto_submit:
        os.system(f'condor_submit {sub_path}')
    else:
        print(f'auto_submit is False. Submit with:\n  condor_submit {sub_path}')
