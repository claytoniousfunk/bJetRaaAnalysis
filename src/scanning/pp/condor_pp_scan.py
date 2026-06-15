import os
import math

jobname = 'pp_HighEGJet_Jet15_scan'

#put the file list here
#dblist = '../../../fileNames/fileNames_pp_SingleMuon_24Oct23.txt'
#dblist = '../../../fileNames/fileNames_pp_MinBias.txt'
dblist = '../../../fileNames/fileNames_pp_HighEGJet.txt'

exe = 'pp_scan.C'
time_flavour = '"workday"'  # 8h
nsplit = 5  # input files per condor job

# -----------------------------------------------------------------------

pwd = os.getenv('PWD')
files = [l.strip() for l in open(dblist).readlines() if l.strip()]
nfiles = len(files)
njobs = int(math.ceil(float(nfiles) / nsplit))

os.makedirs(jobname, exist_ok=True)
os.makedirs(f'{jobname}/log', exist_ok=True)

sub_blocks = [
    'Universe = vanilla',
    f'Log        = {pwd}/{jobname}/log/job_$(ProcId).log',
    f'Output     = {pwd}/{jobname}/log/job_$(ProcId).out',
    f'Error      = {pwd}/{jobname}/log/job_$(ProcId).err',
    'getenv     = True',
    'x509userproxy = $ENV(X509_USER_PROXY)',
    'use_x509userproxy = True',
    f'+JobFlavour = {time_flavour}',
    '',
]

for i in range(njobs):
    start = i * nsplit + 1        # 1-indexed group numbers for pp_scan.C
    end   = min((i + 1) * nsplit, nfiles)

    script = f'#!/bin/bash\ncd {pwd}\n'
    for idx in range(start, end + 1):
        script += f"root -l -b -q '{exe}({idx})'\n"

    script_path = f'{pwd}/{jobname}/script_{i}.sh'
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)

    sub_blocks.append(f'Executable = {script_path}')
    sub_blocks.append('Queue')
    sub_blocks.append('')

sub_path = f'{jobname}/condor_submit.cfg'
with open(sub_path, 'w') as f:
    f.write('\n'.join(sub_blocks))

print(f'Generated {njobs} jobs for {nfiles} input files (nsplit={nsplit})')
os.system(f'condor_submit {sub_path}')
