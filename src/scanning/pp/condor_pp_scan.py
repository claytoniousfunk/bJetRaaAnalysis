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

# Scripts and logs go to AFS home so the schedd can reach them
# regardless of whether the repo is on EOS or AFS.
afs_job_dir = os.path.join(os.getenv('HOME'), 'condor_jobs', jobname)

# -----------------------------------------------------------------------

src_dir = os.getenv('PWD')  # where pp_scan.C lives (may be on EOS)
files = [l.strip() for l in open(dblist).readlines() if l.strip()]
nfiles = len(files)
njobs = int(math.ceil(float(nfiles) / nsplit))

os.makedirs(afs_job_dir, exist_ok=True)
os.makedirs(os.path.join(afs_job_dir, 'log'), exist_ok=True)

# generate one script per job
for i in range(njobs):
    start = i * nsplit + 1        # 1-indexed group numbers for pp_scan.C
    end   = min((i + 1) * nsplit, nfiles)

    script = f'#!/bin/bash\ncd {src_dir}\n'
    for idx in range(start, end + 1):
        script += f"root -l -b -q '{exe}({idx})'\n"

    script_path = os.path.join(afs_job_dir, f'script_{i}.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)

# single submit file with one Queue statement
sub = f"""Universe   = vanilla
Executable = {afs_job_dir}/script_$(ProcId).sh
Log        = {afs_job_dir}/log/job_$(ProcId).log
Output     = {afs_job_dir}/log/job_$(ProcId).out
Error      = {afs_job_dir}/log/job_$(ProcId).err
getenv     = True
x509userproxy = $ENV(X509_USER_PROXY)
use_x509userproxy = True
+JobFlavour = {time_flavour}
Queue {njobs}
"""

sub_path = os.path.join(afs_job_dir, 'condor_submit.cfg')
with open(sub_path, 'w') as f:
    f.write(sub)

print(f'Generated {njobs} jobs for {nfiles} input files (nsplit={nsplit})')
print(f'Scripts/logs → {afs_job_dir}')
os.system(f'condor_submit {sub_path}')
