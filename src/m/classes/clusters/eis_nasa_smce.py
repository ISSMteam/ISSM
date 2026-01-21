import os
import shutil
import subprocess

try:
    from eis_nasa_smce_settings import eis_nasa_smce_settings
except ImportError:
    print('You need eis_nasa_smce_settings.py to proceed, check presence and sys.path')
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmssh import issmssh
from MatlabFuncs import *
from pairoptions import pairoptions

class eis_nasa_smce(object):
    """EIS_NASA_SMCE cluster class definition

    Usage:
        cluster = eis_nasa_smce()
        cluster = eis_nasa_smce('np', 3)
        cluster = eis_nasa_scme('np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = '52.10.233.96'
        self.login = 'jdquinn1'
        self.idfile = '~/.ssh/eis-nasa-smce'
        self.modules = ['intelmpi']
        self.numnodes = 4
        self.cpuspernode = 1
        self.port = 0
        self.time = 12 * 60 * 60
        self.processor = 'skylake'
        self.partition = 'sealevel-c5xl-spot'
        self.srcpath = '/efs/issm-new/binaries/repos/trunk-jpl-working'
        self.extpkgpath = '/efs/issm-new/binaries/ext'
        self.codepath = '/efs/issm-new/binaries/repos/trunk-jpl-working/bin'
        self.executionpath = '~/issm-exec'
        self.interactive = 0
        self.numstreams = 1
        self.hyperthreading = 0
        self.email = ''

        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        try:
            self = eis_nasa_smce_settings(self)
        except NameError:
            print('eis_nasa_smce_settings.py not found, using default settings')

        # OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class eis_nasa_smce object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    idfile: {}\n'.format(self.idfile)
        s += '    modules: {}\n'.format(strjoin(self.modules, ', '))
        s += '    numnodes: {}\n'.format(self.numnodes)
        s += '    cpuspernode: {}\n'.format(self.cpuspernode)
        s += '    np: {}\n'.format(self.nprocs())
        s += '    port: {}\n'.format(self.port)
        s += '    time: {}\n'.format(self.time)
        s += '    processor: {}\n'.format(self.processor)
        s += '    partition: {}\n'.format(self.partition)
        s += '    srcpath: {}\n'.format(self.srcpath)
        s += '    extpkgpath: {}\n'.format(self.extpkgpath)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    interactive: {}\n'.format(self.interactive)
        s += '    numstreams: {}\n'.format(self.numstreams)
        s += '    hyperthreading: {}\n'.format(self.hyperthreading)
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.numnodes * self.cpuspernode
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Now, check cluster.cpuspernode according to processor type
        #if self.processor == 'skylake':
        #    if self.cpuspernode > 14 or self.cpuspernode < 1:
        #        md = md.checkmessage('cpuspernode should be between 1 and 14 for \'skyw\' processors in hyperthreading mode')
        #else:
        #    md = md.checkmessage('unknown processor type, should be \'skylake\'')

        # Miscellaneous
        if not self.login:
            md = md.checkmessage('login empty')
        if self.port:
            md = md.checkmessage('port must be set to 0 as we do not have an SSH tunnel')
        if not self.codepath:
            md = md.checkmessage('codepath empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath empty')

        return self
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        if isgprof:
            print('gprof not supported by cluster, ignoring...')

        issmexec = 'issm.exe'
        mpiexec = 'mpiexec' # Set to alternative mpiexec if desired

        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')[0:2]
            version = float(str(version[0]))
            if version >= 6:
                issmexec = 'issm_dakota.exe'
        if isoceancoupling:
            issmexec = 'issm_ocean.exe'

        # Write queuing script
        fid = open(modelname + '.queue', 'w')

        fid.write('#!/bin/bash\n')
        fid.write('#SBATCH --partition={} \n'.format(self.partition))
        fid.write('#SBATCH -J {} \n'.format(modelname))
        fid.write('#SBATCH -o {}.outlog \n'.format(modelname))
        fid.write('#SBATCH -e {}.errlog \n'.format(modelname))
        fid.write('#SBATCH --nodes={} \n'.format(self.numnodes))
        fid.write('#SBATCH --ntasks-per-node={} \n'.format(self.cpuspernode))
        fid.write('#SBATCH --cpus-per-task={} \n'.format(self.numstreams))
        fid.write('#SBATCH -t {:02d}:{:02d}:00 \n'.format(int(floor(self.time / 3600)), int(floor(self.time % 3600) / 60)))
        if (self.email.find('@')>-1):
            fid.write('#SBATCH --mail-user={} \n'.format(self.email))
            fid.write('#SBATCH --mail-type=BEGIN,END,FAIL \n\n')
        fid.write('source /etc/profile\n')
        fid.write('source /shared/spack/share/spack/setup-env.sh\n')
        for i in range(len(self.modules)):
             fid.write('module load {} &> /dev/null\n'.format(self.modules[i]))
        fid.write('export MPI_GROUP_MAX=64\n\n')
        fid.write('export MPI_UNBUFFERED_STDIO=true\n\n')
        fid.write('export PATH="$PATH:/opt/slurm/bin"\n')
        fid.write('export PATH="$PATH:."\n\n')
        fid.write('export ISSM_DIR="{}"\n'.format(self.srcpath))
        fid.write('export ISSM_EXT_DIR="{}"\n'.format(self.extpkgpath))
        fid.write('source $ISSM_DIR/etc/environment.sh\n')
        fid.write('cd {}/{}/\n\n'.format(self.executionpath, dirname))
        fid.write('{} -n {} {}/{} {} {}/{} {}\n'.format(mpiexec, self.nprocs(), self.codepath, issmexec, solution, self.executionpath, dirname, modelname))

        if not io_gather: # concatenate the output files
            fid.write('cat {}.outbin.* > {}.outbin'.format(modelname, modelname))
        fid.close()

        # In interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(modelname + '.run', 'w')
            if not isvalgrind:
                fid.write('{} -np {} {}/{} {} {}/{} {}\n'.format(mpiexec, self.nprocs(), self.codepath, issmexec, solution, self.executionpath, dirname, modelname))
            else:
                fid.write('{} -np {} valgrind --leak-check=full {}/{} {} {}/{} {}\n'.format(mpiexec, self.nprocs(), self.codepath, issmexec, solution, self.executionpath, dirname, modelname))
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

        print('uploading input file and queueing script')
        if self.interactive:
            directory = '{}/Interactive{}'.format(self.executionpath, self.interactive)
        else:
            directory = self.executionpath

        # NOTE: Replacement for issmscpout(self.name, directory, self.login, self.port, ['{}.tar.gz'.format(dirname)])
        copystring = 'cp {}.tar.gz /efs/issm/tmp'.format(dirname)
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
                launchcommand = 'cd {} && rm -rf {} && mkdir {} && cd {} && cp /efs/issm/tmp/{}.tar.gz . && tar -zxf {}.tar.gz && /opt/slurm/bin/sbatch {}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)

        #Execute Queue job

        # NOTE: Replacement for issmssh(self.name, self.login, self.port, launchcommand)
        subprocess.call('ssh -l {} -i {} {} "{}"'.format(self.login, self.idfile, self.name, launchcommand), shell=True)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
    
        # NOTE: Replacement for issmscpin(self.name, self.login, self.port, directory, filelist)
        directory = '{}/{}/'.format(self.executionpath, dirname)
        fileliststr = '{' + ','.join([str(x) for x in filelist]) + '}'
        downloadcommand = 'scp -i {} {}@{}:{} {}/.'.format(self.idfile, self.login, self.name, os.path.join(directory, fileliststr), os.getcwd())
        subprocess.call(downloadcommand, shell=True) 
    # }}}
