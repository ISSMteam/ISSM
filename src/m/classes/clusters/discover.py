import subprocess

try:
    from discover_settings import discover_settings
except ImportError:
    print('You need discover_settings.py to proceed, check presence and sys.path')
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from MatlabFuncs import *
from pairoptions import pairoptions
from QueueRequirements import QueueRequirements

class discover(object):
    """DISCOVER cluster class definition

    Usage:
        cluster = discover()
        cluster = discover('np', 3)
        cluster = discover('np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = 'login.nccs.nasa.gov'
        self.login = ''
        self.modules = ['comp/intel/20.0.0.166', 'mpi/sgi-mpt/2.17', 'cmake/3.17.0']
        self.numnodes = 20
        self.cpuspernode = 8
        self.port = 0
        self.queue = 'general'
        self.time = 12 * 60 * 60
        self.processor = 'west'
        self.srcpath = ''
        self.codepath = ''
        self.executionpath = ''
        self.grouplist = ''
        self.interactive = 0
        self.bbftp = 0
        self.numstreams = 8
        self.hyperthreading = 0
        self.email = ''

        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        try:
            self = discover_settings(self)
        except NameError:
            print('discover_settings.py not found, using default settings')

        # OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class discover object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    modules: {}\n'.format(strjoin(self.modules, ', '))
        s += '    numnodes: {}\n'.format(self.numnodes)
        s += '    cpuspernode: {}\n'.format(self.cpuspernode)
        s += '    np: {}\n'.format(self.nprocs())
        s += '    port: {}\n'.format(self.port)
        s += '    queue: {}\n'.format(self.queue)
        s += '    time: {}\n'.format(self.time)
        s += '    processor: {}\n'.format(self.processor)
        s += '    srcpath: {}\n'.format(self.srcpath)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    grouplist: {}\n'.format(self.grouplist)
        s += '    interactive: {}\n'.format(self.interactive)
        s += '    bbftp: {}\n'.format(self.bbftp)
        s += '    numstreams: {}\n'.format(self.numstreams)
        s += '    hyperthreading: {}\n'.format(self.hyperthreading)
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.numnodes * self.cpuspernode
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        queuedict = {'long': [24 * 60 * 60, 560],
                     'allnccs': [12 * 60 * 60, 6000],
                     'debug': [1 * 60 * 60, 532]}
        QueueRequirements(queuedict, self.queue, self.nprocs(), self.time)

        # Now, check cluster.cpuspernode according to processor type
        if self.processor == 'sand':
            if self.cpuspernode > 16 or self.cpuspernode < 1:
                md = md.checkmessage('cpuspernode should be between 1 and 16 for \'sand\' processors in hyperthreading mode')
        elif self.processor == 'hasw':
            if self.cpuspernode > 28 or self.cpuspernode < 1:
                md = md.checkmessage('cpuspernode should be between 1 and 28 for \'hasw\' processors in hyperthreading mode')
        else:
            md = md.checkmessage('unknown processor type, should be \'sand\' or \'hasw\'')

        # Miscellaneous
        if not self.login:
            md = md.checkmessage('login empty')
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

        fid.write('#!/bin/bash\n')
        fid.write('#SBATCH -J {} \n'.format(modelname))
        fid.write('#SBATCH --qos={} \n'.format(self.queue))
        fid.write('#SBATCH -o {}.outlog \n'.format(modelname))
        fid.write('#SBATCH -e {}.errlog \n'.format(modelname))
        fid.write('#SBATCH -n {} \n'.format(self.nprocs()))
        fid.write('#SBATCH -N {} \n'.format(self.numnodes))
        fid.write('#SBATCH -t {:02d}:{:02d}:00 \n'.format(int(floor(self.time / 3600)), int(floor(self.time % 3600) / 60)))
        fid.write('#SBATCH -A {} \n\n'.format(self.grouplist))
        if (self.email.find('@')>-1):
            fid.write('#SBATCH --mail-user={} \n'.format(self.email))
            fid.write('#SBATCH --mail-type=end \n\n')
        fid.write('. /usr/share/modules/init/bash\n\n')
        for i in range(len(self.modules)):
            fid.write('module load {}\n'.format(self.modules[i]))
        fid.write('export MPI_GROUP_MAX=64\n\n')
        fid.write('export MPI_UNBUFFERED_STDIO=true\n\n')
        fid.write('export PATH="$PATH:."\n\n')
        fid.write('export ISSM_DIR="{}/../"\n'.format(self.codepath)) # FIXME
        fid.write('source $ISSM_DIR/etc/environment.sh\n') # FIXME
        fid.write('cd {}/{}\n\n'.format(self.executionpath, dirname))

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

        print('uploading input file and queuing script')
        if self.interactive:
            directory = '{}/Interactive{}'.format(self.executionpath, self.interactive)
        else:
            directory = self.executionpath

        if self.bbftp:
            issmbbftpout(self.name, directory, self.login, self.port, self.numstreams, '{}.tar.gz'.format(dirname))
        else:
            issmscpout(self.name, directory, self.login, self.port, ['{}.tar.gz'.format(dirname)])
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

        print('launching solution sequence on remote cluster')
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        if self.interactive:
            directory = '{}/Interactive{}'.format(self.executionpath, self.interactive)
        else:
            directory = '{}/{}/'.format(self.executionpath, dirname)

        if self.bbftp:
            issmbbftpin(self.name, self.login, self.port, self.numstreams, directory, filelist)
        else:
            issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
