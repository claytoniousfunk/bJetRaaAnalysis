import math
import subprocess
import os

eos_purdue_t2 = 'root://xrootd.rcac.purdue.edu/'
eos_dir_list = {'cern': '/eos/cms/store/group/phys_heavyions/wangx',
                'purdue': '', '': '', 'lpc': 'root://cmseos.fnal.gov//eos/uscms/store/user/wangx'}
time_sequence = {'none': '', '20m': '"espresso"', '1h': '"microcentury"', '2h': '"longlunch"',
                 '8h': '"workday"', '1d': '"tomorrow"', '3d': '"testmatch"', '1w': '"nextweek"'}

_here = os.path.dirname(os.path.abspath(__file__))


def makeTdrBall(eospath, make=False):
    eospath = eospath + 'tarball/'
    cmssw_dir = os.getenv('CMSSW_BASE')
    cmssw_ver = os.getenv('CMSSW_VERSION')
    if make:
        subprocess.call(['tar', '-zvcf', cmssw_ver + '.tgz', '-C', cmssw_dir, '.'])
        subprocess.call(['xrdcp', '-p', '-f', cmssw_ver + '.tgz', eospath + cmssw_ver + '.tgz'])
    return ('mkdir ' + cmssw_ver + '\npushd ' + cmssw_ver + '\nxrdcp -f ' + eospath + cmssw_ver +
            '.tgz ./\ntar -xf ' + cmssw_ver + '.tgz\npushd ' + cmssw_ver + '/src\nscramv1 b ProjectRename\n')


def subText(inf, outf, parlist):
    if len(parlist) == 0:
        print('Error: Empty list abort!')
        return
    ff = open(inf, 'r')
    lines = ff.readlines()
    of = open(outf, 'w+')
    newlines = []
    for line in lines:
        newline = line
        for k in parlist.keys():
            newline = newline.replace(k, str(parlist[k]), 1)
        newlines.append(newline)
    of.writelines(newlines)
    of.close()


def mkdirCheck(workfolder):
    if not os.path.exists('{FOLDER}'.format(FOLDER=workfolder)):
        os.system('mkdir -p -v {FOLDER}'.format(FOLDER=workfolder))
    else:
        print('FD exits! Abord!'.format(FD=workfolder))
        return


class jobManager:
    def __init__(self, jobSite, jobname, method, executable, runlist, output_dir='', infile=[],
                 env_mode='local', time='20m', remakeTarball=False, rm_list=[]):
        self.rm_list = rm_list
        self.jobname = jobname
        self.jobSite = jobSite
        self.executable = executable
        self.method = method
        self.runlist = runlist
        self.output_dir = output_dir
        self.cfg_template = 'default'
        self.cfg_list = []
        self.env_mode = env_mode
        self.store_path = eos_dir_list[jobSite]
        self.time = time
        self.make_tarball = remakeTarball
        self.inputFiles = infile
        self.nsplit = 1
        self.binary = ''

    def set_env_prefix(self):
        env = 'export PYTHONHOME="/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15/"\n'
        return env

    def set_run_environment_local(self):
        pwd = os.getenv('PWD')
        base = os.getenv('CMSSW_BASE')
        if base is None:
            # no CMSSW environment — plain ROOT job, just cd to the job directory
            cmd = 'pushd ' + pwd + '/' + self.jobname + '\n'
        else:
            cmd = 'pushd ' + base + '/src\n' + 'eval `scramv1 runtime -csh`\npushd '
            cmd = cmd + pwd + '/' + self.jobname + '\n'
        return cmd

    def set_run_environment_tarball(self):
        path = self.store_path
        cmd = makeTdrBall(path, self.make_tarball)
        pwd = os.getenv('PWD')
        ver = os.getenv('CMSSW_VERSION')
        cmd = cmd + 'eval `scramv1 runtime -csh`\n'
        working_dir = pwd.split(ver + '/src')
        cmd = cmd + 'popd\n'
        cmd = cmd + 'pushd .' + working_dir[1]
        return cmd

    def set_run_environment(self):
        if self.env_mode == 'tarball':
            return self.set_run_environment_tarball()
        elif self.env_mode == 'local':
            return self.set_run_environment_local()

    def set_run_time(self):
        time_s = '+JobFlavour = ' + time_sequence[self.time] + '\n'
        return time_s

    def generate_cfg(self):
        nsplit = self.nsplit
        workfolder = self.jobname
        if not os.path.exists('{FOLDER}'.format(FOLDER=workfolder)):
            os.system('mkdir {FOLDER}'.format(FOLDER=workfolder))
        os.system('mkdir {FOLDER}/outCondor'.format(FOLDER=workfolder))
        os.system('mkdir {FOLDER}/data'.format(FOLDER=workfolder))
        os.system('cp -v {exe} {FOLDER}/'.format(FOLDER=workfolder, exe=self.executable))
        for files in self.inputFiles:
            os.system('cp -v {exe} {FOLDER}/'.format(FOLDER=workfolder, exe=files))
        if self.method == 'root':
            self.binary = 'root -b -l -q'
        elif self.method == 'cmsRun':
            self.binary = 'cmsRun'
        self.cmsswDir = os.getenv('CMSSW_BASE')
        pwd = os.getenv('PWD')
        files = open(self.runlist).readlines()
        outputname0 = 'job_output'
        file_keep = ''
        njobs = int(math.ceil(float(len(files)) / nsplit))
        prerun_cmd = ''
        env_setup = self.set_env_prefix()
        env_setup = env_setup + self.set_run_environment()
        output_cmd = 'cp -v data/' + file_keep + ' $LS_SUBCWD' + '/data/' + file_keep
        if self.env_mode == 'local':
            output_cmd = ''

        # use templates bundled alongside this file
        script_temp = os.path.join(_here, 'utility', 'script_condor_template.sh')
        cfg_temp = os.path.join(_here, 'utility', 'condor_template.cfg')

        ifiles = 0
        outputname = ''
        for i in range(njobs):
            if nsplit != 1:
                outputname = './temp' + str(i) + '/' + outputname0
                prerun_cmd = 'mkdir temp' + str(i) + '\n'
            else:
                outputname = './data/' + outputname0
            starti = i * nsplit
            endi = (i + 1) * nsplit
            cmdline = ''
            for f in files[starti:endi]:
                if self.method == 'root':
                    cmdline = (cmdline + self.binary + ' ' + self.executable +
                               '\\("' + f.rstrip() + '","' + outputname + '_' + str(ifiles) + '.root")\\\n')
                else:
                    cmdline = (cmdline + self.binary + ' ' + self.executable + ' ' +
                               f.rstrip() + ' ' + outputname + '_' + str(ifiles) + '.root\n')
                for frm in self.rm_list:
                    cmdline = cmdline + 'rm ' + frm.format(job_number=ifiles) + '\n'
                ifiles = ifiles + 1
            if nsplit != 1:
                file_keep = 'job_output_part' + str(i) + '.root'
                cmdline = cmdline + 'hadd -f data/' + file_keep + ' ./temp' + str(i) + '/*.root\n'
                cmdline = cmdline + 'rm -r temp' + str(i) + '\n'
            else:
                file_keep = 'data/' + outputname + str(i) + '.root'
            parlist = {'EXECUTABLE': cmdline}
            parlist['ENV_SETUP'] = env_setup
            parlist['PRERUN'] = prerun_cmd
            parlist['OUTPUT'] = output_cmd
            subText(script_temp, workfolder + '/script_' + str(i) + '.sh', parlist)
            parlist = {}
            parlist['SCRIPT'] = 'script_' + str(i) + '.sh'
            os.system('chmod 755 ' + workfolder + '/script_' + str(i) + '.sh')
            parlist['KINDEX'] = str(i)
            parlist['FOLDER'] = pwd + '/' + workfolder
            parlist['TIME'] = self.set_run_time()
            subText(cfg_temp, workfolder + '/condor_cfg_' + str(i) + '.cfg', parlist)
            self.cfg_list.append('condor_cfg_' + str(i) + '.cfg')

    def submit(self):
        os.chdir(self.jobname)
        for cfg in self.cfg_list:
            os.system("condor_submit {CFG}".format(CFG=cfg))
