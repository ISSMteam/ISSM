from datetime import datetime
import os
from random import randint
import subprocess

try:
    from ub_ccr_settings import ub_ccr_settings
except ImportError:
    print('You need ub_ccr_settings.py to proceed, check presence and sys.path')
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from MatlabFuncs import *
from pairoptions import pairoptions
from QueueRequirements import QueueRequirements

class ub_ccr(object):
    """ub_ccr cluster class definition

    Usage:
        cluster = ub_ccr()
        cluster = ub_ccr('np', 3)
        cluster = ub_ccr('np', 3, 'login', 'username')

    See Also: https://docs.ccr.buffalo.edu/en/latest/hpc/jobs/#requesting-cores-and-nodes
    """

    def __init__(self, *args):  # {{{
        self.name = oshostname()
        self.login = ''
        self.port = ''
        self.cluster = 'ub-hpc'
        self.partition = 'general-compute'
        self.qos = 'general-compute'
        self.account = ''
        self.time = 1 * 60 * 60
        self.ntasks = 1
        self.cpuspertask = 1
        self.exclusive = False
        self.mem = '10000'
        self.jobname = ''
        self.modules = ['ccrsoft/2023.01', 'gcc/11.2.0']
        self.srcpath = '/projects/grid/ghub/ISSM/repos/ISSM'
        self.codepath = '/projects/grid/ghub/ISSM/repos/ISSM/bin'
        self.executionpath = '$HOME/ISSM/execution'
        self.interactive = 0
        self.bbftp = 0
        self.email = ''

        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        try:
            self = ub_ccr_settings(self)
        except NameError:
            print('ub_ccr_settings.py not found, using default settings')

        # Get other fields
        self = options.AssignObjectFields(self)

        # Set jobname if empty
        now = datetime.now()
        now_str = now.strftime('%Y%m%d%H%M%S')
        rand_int = randint(1000, 9999)
        jobname = '{}-issm-{}-{}'.format(self.login, now_str, rand_int)
        self.jobname = jobname
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class ub_ccr object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    port: {}\n'.format(self.port)
        s += '    cluster: {}\n'.format(self.cluster)
        s += '    partition: {}\n'.format(self.partition)
        s += '    qos: {}\n'.format(self.qos)
        s += '    account: {}\n'.format(self.account)
        s += '    time: {}\n'.format(self.time)
        s += '    ntasks: {}\n'.format(self.ntasks)
        s += '    cpuspertask: {}\n'.format(self.cpuspertask)
        s += '    exclusive: {}\n'.format(self.exclusive)
        s += '    mem: {}\n'.format(self.mem)
        s += '    jobname: {}\n'.format(self.jobname)
        s += '    modules: {}\n'.format(strjoin(self.modules, ', '))
        s += '    srcpath: {}\n'.format(self.srcpath)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    grouplist: {}\n'.format(self.grouplist)
        s += '    interactive: {}\n'.format(self.interactive)
        s += '    bbftp: {}\n'.format(self.bbftp)
        s += '    login: {}\n'.format(self.email)
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.ntasks * self.cpuspertask
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if self.cluster == 'ub-hpc':
            if self.qos != self.partition:
                if self.qos not in ['supporters', 'mri', 'nih']:
                    md = md.checkmessage('Value of \'qos\' should either match value of \'partition\' or be set to one of [\'supporters\', \'mri\', \'nih\']')
            queuedict = {
                'debug': [1 * 60 * 60, 64],
                'general-compute': [72 * 60 * 60, 64],
                'industry': [6 * 60 * 60, 56],
                'scavenger': [6 * 60 * 60, 64],
                'viz': [24 * 60 * 60, 64],
                'faculty': [72 * 60 * 60, 64],
            }
            QueueRequirements(queuedict, self.partition, self.nprocs(), self.time)
        elif self.cluster == 'faculty':
            # TODO: Add checks for max values based on particular partition
            if self.account == '':
                md = md.checkmessage('please supply valid \'account\' when using \'faculty\' cluster')
            if self.partition != 'sophien' and self.qos != 'sophien' and self.account != 'sophien':
                md = md.checkmessage('combination of \'partition\' and \'qos\' and \'account\' invalid')
        else:
            md = md.checkmessage('invalid value for \'cluster\'')

        # Miscellaneous
        if not self.srcpath:
            md = md.checkmessage('srcpath empty')
        if not self.codepath:
            md = md.checkmessage('codepath empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath empty')

        return self
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        if isgprof:
            print('gprof not supported by cluster, ignoring...')

        executable = 'issm.exe'
        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')[0:2]
            version = float(str(version[0]))
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Write queuing script
        fid = open(modelname + '.queue', 'w')

        fid.write('#!/bin/bash -l\n')
        fid.write('#SBATCH --time {:02d}:{:02d}:00\n'.format(int(floor(self.time / 3600)), int(floor(self.time % 3600) / 60)))
        fid.write('#SBATCH --ntasks={}\n'.format(self.ntasks))
        fid.write('#SBATCH --cpus-per-task={}\n'.format(self.cpuspertask))
        if self.exclusive:
            fid.write('#SBATCH --exclusive\n')
        if (self.cluster != 'faculty'):
            fid.write('#SBATCH --constraint="[SAPPHIRE-RAPIDS-IB|ICE-LAKE-IB|CASCADE-LAKE-IB|EMERALD-RAPIDS-IB]"\n')
        fid.write('#SBATCH --mem={}\n'.format(self.mem))
        fid.write('#SBATCH --job-name={}\n'.format(self.jobname))
        fid.write('#SBATCH --output={}.outlog\n'.format(modelname))
        fid.write('#SBATCH --error={}.errlog\n'.format(modelname))
        if (self.email.find('@')>-1):
            fid.write('#SBATCH --mail-user={}\n'.format(self.email))
            fid.write('#SBATCH --mail-type=end\n')
        fid.write('#SBATCH --partition={}\n'.format(self.partition))
        fid.write('#SBATCH --qos={}\n'.format(self.qos))
        if (self.account != ''):
            fid.write('#SBATCH --account={}\n'.format(self.account))
        fid.write('#SBATCH --cluster={}\n'.format(self.cluster))

        #fid.write('. /usr/share/modules/init/bash\n\n')
        for i in range(len(self.modules)):
            fid.write('module load {}\n'.format(self.modules[i]))
        fid.write('export PATH="$PATH:."\n\n')
        fid.write('export ISSM_DIR="{}"\n'.format(self.srcpath))
        fid.write('source $ISSM_DIR/etc/environment.sh\n')
        fid.write('cd {}/{}\n\n'.format(self.executionpath, dirname))

        fid.write('which mpiexec\n\n')

        fid.write('mpiexec -np {} {}/{} {} {}/{} {}\n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, dirname, modelname))

        if not io_gather: # concatenate the output files
            fid.write('cat {}.outbin.* > {}.outbin'.format(modelname, modelname))
        fid.close()

        # In interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(modelname + '.run', 'w')
            if not isvalgrind:
                fid.write('mpiexec -np {} {}/{} {} {}/{} {}\n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, dirname, modelname))
            else:
                fid.write('mpiexec -np {} valgrind --leak-check=full {}/{} {} {}/{} {}\n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, dirname, modelname))
            if not io_gather: # concatenate the output files
                fid.write('cat {}.outbin.* > {}.outbin'.format(modelname, modelname))
            fid.close()
            fid = open(modelname + '.errlog', 'w') # TODO: Change this to system call (touch <file>)?
            fid.close()
            fid = open(modelname + '.outlog', 'w') # TODO: Change this to system call (touch <file>)?
            fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one zip
        compressstring = 'tar -zcf {}.tar.gz'.format(dirname)
        for file in filelist:
            compressstring += ' {}'.format(file)
        if self.interactive:
            compressstring += ' {}.run {}.errlog {}.outlog'.format(modelname, modelname, modelname)
        subprocess.call(compressstring, shell=True)

        #upload input files
        if self.interactive:
            directory = '{}/Interactive{}'.format(self.executionpath, self.interactive)
        else:
            directory = self.executionpath

        # if self.bbftp:
        #     issmbbftpout(self.name, directory, self.login, self.port, self.numstreams, '{}.tar.gz'.format(dirname))
        # else:
        #     issmscpout(self.name, directory, self.login, self.port, ['{}.tar.gz'.format(dirname)])

        # NOTE: Replacement for issmscpout(self.name, directory, self.login, self.port, ['{}.tar.gz'.format(dirname)])
        mkdirstring = 'if [ ! -d {} ]; then mkdir -p {}; fi'.format(self.executionpath, self.executionpath)
        subprocess.call(mkdirstring, shell=True)
        copystring = 'cp {}.tar.gz {}'.format(dirname, self.executionpath)
        subprocess.call(copystring, shell=True)
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        if self.interactive:
            if not isempty(restart):
                launchcommand = 'cd {}/Interactive{}'.format(self.executionpath, self.interactive)
            else:
                launchcommand = 'cd {}/Interactive{} && tar -zxf {}.tar.gz'.format(self.executionpath, self.interactive, dirname)
        else:
            if not isempty(restart):
                launchcommand = 'cd {} && cd {} && sbatch {}.queue'.format(self.executionpath, dirname, modelname)
            else:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && sbatch {}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)

        #Execute Queue job
        # NOTE: Replacement for issmssh(self.name, self.login, self.port, launchcommand)
        subprocess.call(launchcommand, shell=True)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        if self.interactive:
            directory = '{}/Interactive{}'.format(self.executionpath, self.interactive)
        else:
            directory = '{}/{}/'.format(self.executionpath, dirname)

        # if self.bbftp:
        #     issmbbftpin(self.name, self.login, self.port, self.numstreams, directory, filelist)
        # else:
        #     issmscpin(self.name, self.login, self.port, directory, filelist)

        # NOTE: Replacement for issmscpin(self.name, self.login, self.port, directory, filelist)
        filelist = [os.path.join(directory, x) for x in filelist]
        fileliststr = ' '.join([str(x) for x in filelist])
        downloadcommand = 'cp {} {}'.format(fileliststr, os.getcwd())
        subprocess.call(downloadcommand, shell=True)
    # }}}
