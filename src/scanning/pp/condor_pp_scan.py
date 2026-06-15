import os
import math

jobname = 'pp_HighEGJet_Jet100_scan'

#put the file list here
#dblist = '../../../fileNames/fileNames_pp_SingleMuon_24Oct23.txt'
#dblist = '../../../fileNames/fileNames_pp_MinBias.txt'
dblist = '../../../fileNames/fileNames_pp_HighEGJet.txt'

exe = 'pp_scan.C'
time_flavour = '"workday"'  # 8h
nsplit = 5  # input files per condor job

# -----------------------------------------------------------------------
# Determine where to put scripts and logs.
#
# The condor standard batch schedds at CERN reject /eos paths in submit
# files. Prefer AFS (which they can always reach). If AFS is unavailable
# (user migrated fully to EOS), fall back to EOS and hope the user has
# enabled the EosSubmit schedd or will move scripts manually.
#
# Run this on lxplus to see which branch is taken:
#   python3 condor_pp_scan.py 2>&1 | head -5
# -----------------------------------------------------------------------

_user = os.getenv('USER') or os.getenv('LOGNAME') or ''
_afs  = f'/afs/cern.ch/user/{_user[0]}/{_user}' if _user else ''

if _afs and os.path.isdir(_afs):
    job_dir = os.path.join(_afs, 'condor_jobs', jobname)
    eos_submit_flag = ''
    print(f'[condor] Using AFS for scripts/logs: {job_dir}')
else:
    job_dir = os.path.join(os.getenv('HOME', f'/tmp/{_user}'), 'condor_jobs', jobname)
    eos_submit_flag = '+EOS = True\n'
    print(f'[condor] AFS not found — using EOS (needs EosSubmit schedd): {job_dir}')

# -----------------------------------------------------------------------

src_dir = os.getenv('PWD')  # repo may be on EOS; that's fine for the scripts
files   = [l.strip() for l in open(dblist).readlines() if l.strip()]
nfiles  = len(files)
njobs   = int(math.ceil(float(nfiles) / nsplit))

os.makedirs(job_dir, exist_ok=True)
os.makedirs(os.path.join(job_dir, 'log'), exist_ok=True)

# generate one script per job
for i in range(njobs):
    start = i * nsplit + 1
    end   = min((i + 1) * nsplit, nfiles)

    script = f'#!/bin/bash\ncd {src_dir}\n'
    for idx in range(start, end + 1):
        script += f"root -l -b -q '{exe}({idx})'\n"

    script_path = os.path.join(job_dir, f'script_{i}.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)

# single submit file with one Queue statement
sub = f"""Universe   = vanilla
Executable = {job_dir}/script_$(ProcId).sh
Log        = {job_dir}/log/job_$(ProcId).log
Output     = {job_dir}/log/job_$(ProcId).out
Error      = {job_dir}/log/job_$(ProcId).err
getenv     = True
x509userproxy = $ENV(X509_USER_PROXY)
use_x509userproxy = True
+JobFlavour = {time_flavour}
{eos_submit_flag}Queue {njobs}
"""

sub_path = os.path.join(job_dir, 'condor_submit.cfg')
with open(sub_path, 'w') as f:
    f.write(sub)

print(f'[condor] {njobs} jobs for {nfiles} files (nsplit={nsplit})')
print(f'[condor] submit file → {sub_path}')
os.system(f'condor_submit {sub_path}')
