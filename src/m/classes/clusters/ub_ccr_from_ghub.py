import subprocess

try:
    from ub_ccr_from_ghub_settings import ub_ccr_from_ghub_settings
except ImportError:
    print('You need ub_ccr_from_ghub_settings.py to proceed, check presence and sys.path')
from datetime import datetime
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from MatlabFuncs import *
from pairoptions import pairoptions
from QueueRequirements import QueueRequirements
from random import randint

class ub_ccr_from_ghub(object):
    """ub_ccr_from_ghub cluster class definition

    Usage:
        cluster = ub_ccr_from_ghub()
    """

    def __init__(self, *args):  # {{{
        self.time = 24 * 60 * 60
        self.nodes = 1
        self.ntasks = 1
        self.ntaskspernode = 1
        self.cpuspertask = 8
        self.mem = '8G'
        self.jobname = ''
        self.output = ''
        self.mailuser = ''

        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        try:
            self = ub_ccr_settings(self)
        except NameError:
            print('ub_ccr_settings.py not found, using default settings')

        # OK get other fields
        self = options.AssignObjectFields(self)

        # Override certain fields regardless of how user might have set them
        now = datetime.now()
        now_str = now.strftime('%Y%m%d%H%M%S')
        rand_int = randint(1000, 9999)
        jobname = 'ghub-issm-{}-{}'.format(now_str, rand_int)
        self.jobname = jobname
        self.output = '{}.out'.format(jobname)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class ub_ccr_from_ghub object\n'
        s += '    time: {}\n'.format(self.time)
        s += '    nodes: {}\n'.format(self.nodes)
        s += '    ntasks: {}\n'.format(self.ntasks)
        s += '    ntaskspernode: {}\n'.format(self.ntaskspernode)
        s += '    cpuspertask: {}\n'.format(self.cpuspertask)
        s += '    mem: {}\n'.format(self.mem)
        s += '    jobname: {}\n'.format(self.jobname)
        s += '    output: {}\n'.format(self.output)
        s += '    mailuser: {}\n'.format(self.mailuser)
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.nodes * self.ntaskspernode * self.cpuspertask
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if not (self.ntaskspernode * self.nodes >= self.ntasks):
            md = md.checkmessage('ntasks must be < ntaskspernode * nodes')
        if self.cpuspertask * self.ntaskspernode > 40:
            md = md.checkmessage('maximum number of cores per node (cpuspertask * ntaskspernode) is 40')
        if self.time > 3 * 24 * 60 * 60:
            md = md.checkmessage('time must be >= 72 hours (in seconds)')
        if not self.mem.isdigit():
            if not self.mem[0:-2].isdigit() or not (self.mem[-1] == 'K' or self.mem[-1] == 'M' or self.mem[-1] == 'G' or self.mem[-1] == 'T'):
                md = md.checkmessage('mem should be an integer or be specified with the suffixes: K,M,G,T (default M)')

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

        partition = 'general-compute'
        qos = 'general-compute'
        cluster = 'ub-hpc'
        reservation = 'ubhpc-future'
        modules = ['ccrsoft/2023.01', 'gcc/11.2.0']

        issmprojectdir = '/projects/grid/ghub/ISSM'
        issmdir = '{}/repos/trunk-jpl'
        issmextdir = '{}/ext'
        issmexecdir = '{}/exec'

        fid.write('#!/bin/bash -l\n')
        fid.write('#SBATCH --time={:02d}:{:02d}:00\n'.format(int(floor(self.time / 3600)), int(floor(self.time % 3600) / 60)))
        fid.write('#SBATCH --nodes={}\n'.format(self.nodes))
        fid.write('#SBATCH --ntasks={}\n'.format(self.ntasks))
        fid.write('#SBATCH --ntasks-per-node={}\n'.format(self.ntaskspernode))
        fid.write('#SBATCH --cpus-per-task={}\n'.format(self.cpuspertask))
        fid.write('#SBATCH --mem={}\n'.format(self.mem))
        fid.write('#SBATCH --job-name={}\n'.format(self.jobname))
        fid.write('#SBATCH --output={}.outlog\n'.format(modelname))
        fid.write('#SBATCH --error={}.errlog\n'.format(modelname))
        if (self.mailuser.find('@')>-1):
            fid.write('#SBATCH --mail-user={} \n'.format(self.mailuser))
            fid.write('#SBATCH --mail-type=end \n\n')
        fid.write('#SBATCH --partition={}\n'.format(partition))
        fid.write('#SBATCH --qos={}\n'.format(qos))
        fid.write('#SBATCH --cluster={}\n'.format(cluster))
        fid.write('#SBATCH --reservation={}\n'.format(reservation))
        for i in range(len(modules)):
            fid.write('module load {}\n'.format(modules[i]))
        fid.write('export MPI_GROUP_MAX=64\n\n')
        fid.write('export MPI_UNBUFFERED_STDIO=true\n\n')
        fid.write('export ISSM_DIR="{}"\n'.format(issmdir))
        fid.write('export ISSM_EXT_DIR="{}"\n'.format(issmextdir))
        fid.write('source $ISSM_DIR/etc/environment.sh\n')
        fid.write('cd {}/{}\n\n'.format(issmexecdir, dirname))

        fid.write('mpiexec -np {} {}/{} {} {}/{} {}\n'.format(self.nprocs(), issmdir, executable, solution, issmexecdir, dirname, modelname))

        if not io_gather: # concatenate the output files
            fid.write('cat {}.outbin.* > {}.outbin'.format(modelname, modelname))
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
        directory = issmexecdir

        issmscpout(self.name, directory, self.login, self.port, ['{}.tar.gz'.format(dirname)])
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && sbatch {}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)

        #Execute Queue job
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        directory = '{}/{}/'.format(self.executionpath, dirname)

        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
