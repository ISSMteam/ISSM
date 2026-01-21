import subprocess

from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from MatlabFuncs import *
from pairoptions import pairoptions
try:
    from pfe_settings import pfe_settings
except ImportError:
    # print('You need pfe_settings.py to proceed, check presence and sys.path')
    pass
from QueueRequirements import QueueRequirements

class pfe(object):
    """pfe cluster class definition

    Usage:
        cluster = pfe()
        cluster = pfe('np', 3)
        cluster = pfe('np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = ''
        self.login = ''
        self.modules = ['comp-intel/2018.3.222', 'mpi-intel/2018.3.222', 'scicon/app-tools']
        self.numnodes = 20
        self.cpuspernode = 8
        self.port = 1025
        self.queue = 'long'
        self.time = 12 * 60
        self.processor = 'ivy'
        self.srcpath = ''
        self.extpkgpath = ''
        self.codepath = ''
        self.executionpath = ''
        self.grouplist = ''
        self.interactive = 0
        self.bbftp = 0
        self.numstreams = 8
        self.hyperthreading = 0

        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        try:
            self = pfe_settings(self)
        except NameError:
            #print('pfe_settings.py not found, using default settings')
            pass

        # OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class pfe object\n'
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
        s += '    extpkgpath: {}\n'.format(self.extpkgpath)
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
        queuedict = {'long': [5 * 24 * 60, 2048],
                     'normal': [8 * 60, 2048],
                     'debug': [2 * 60, 150],
                     'devel': [2 * 60, 150]}
        QueueRequirements(queuedict, self.queue, self.nprocs(), self.time)

        # Now, check cluster.cpuspernode according to processor type
        if self.processor == 'ivy':
            if self.hyperthreading:
                if self.cpuspernode > 40 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 40 for \'ivy\' processors in hyperthreading mode')
            else:
                if self.cpuspernode > 20 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 20 for \'ivy\' processors')
        elif self.processor == 'bro':
            if self.hyperthreading:
                if self.cpuspernode > 56 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 56 for \'bro\' processors in hyperthreading mode')
            else:
                if self.cpuspernode > 28 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 28 for \'bro\' processors')
        elif self.processor == 'has':
            if self.hyperthreading:
                if self.cpuspernode > 48 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 48 for \'has\' processors in hyperthreading mode')
            else:
                if self.cpuspernode > 24 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 24 for \'has\' processors')
        elif self.processor == 'san':
            if self.hyperthreading:
                if self.cpuspernode > 32 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 32 for \'san\' processors in hyperthreading mode')
            else:
                if self.cpuspernode > 16 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 16 for \'san\' processors')
        elif self.processor == 'cas_ait':
            if self.hyperthreading:
                if self.cpuspernode > 80 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 80 for \'cas_ait\' processors in hyperthreading mode')
            else:
                if self.cpuspernode > 40 or self.cpuspernode < 1:
                    md = md.checkmessage('cpuspernode should be between 1 and 40 for \'cas_ait\' processors')
        else:
            md = md.checkmessage('unknown processor type, should be \'bro\', \'has\', \'ivy\', \'san\', or \'cas_ait\'')

        # Miscellaneous
        if not self.login:
            md = md.checkmessage('login empty')
        if not self.srcpath:
            md = md.checkmessage('srcpath empty')
        if not self.codepath:
            md = md.checkmessage('codepath empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath empty')
        if not self.grouplist:
            md = md.checkmessage('grouplist empty')

        return self
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        if isgprof:
            print('gprof not supported by cluster, ignoring...')

        executable = 'issm.exe'
        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')[0:2]
            version = float(version)
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Write queuing script
        fid = open(modelname + '.queue', 'w')
        fid.write('#PBS -S /bin/bash\n')
        fid.write('#PBS -l select={}:ncpus={}:model={}\n'.format(self.numnodes, self.cpuspernode, self.processor))
        fid.write('#PBS -l walltime={}\n'.format(self.time * 60)) # walltime is in seconds
        fid.write('#PBS -q {} \n'.format(self.queue))
        fid.write('#PBS -W group_list={}\n'.format(self.grouplist))
        fid.write('#PBS -m e\n')
        fid.write('#PBS -o {}/{}/{}.outlog \n'.format(self.executionpath, dirname, modelname))
        fid.write('#PBS -e {}/{}/{}.errlog \n\n'.format(self.executionpath, dirname, modelname))
        fid.write('. /usr/share/modules/init/bash\n\n')
        for i in range(len(self.modules)):
            fid.write('module load {}\n'.format(self.modules[i]))
        fid.write('export PATH="$PATH:."\n\n')
        fid.write('export MPI_LAUNCH_TIMEOUT=520\n')
        fid.write('export MPI_GROUP_MAX=64\n\n')
        fid.write('export ISSM_DIR="{}"\n'.format(self.srcpath)) # FIXME
        if self.extpkgpath:
            fid.write('export ISSM_EXT_DIR="{}"\n'.format(self.extpkgpath)) 
        fid.write('source $ISSM_DIR/etc/environment.sh\n') # FIXME
        fid.write('cd {}/{}/\n\n'.format(self.executionpath, dirname))
        fid.write('mpiexec -np {} {}/{} {} {}/{} {}\n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, dirname, modelname))

        fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one zip
        compressstring = 'tar -zcf {}.tar.gz'.format(dirname)
        for file in filelist:
            compressstring += ' {}'.format(file)
        subprocess.call(compressstring, shell=True)

        print('uploading input file and queueing script')
        directory = self.executionpath
        issmscpout(self.name, directory, self.login, self.port, [dirname + '.tar.gz'])

    # }}}
    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        # Launch command, to be executed via ssh
        if self.interactive:
            if not isempty(restart):
                launchcommand = 'cd {} /Interactive{}'.format(self.executionpath, self.interactive)
            else:
                if self.interactive == 10:
                    launchcommand = 'cd {}/run && tar -zxf {}.tar.gz'.format(pwd(), dirname)
                else:
                    launchcommand = 'cd {} /Interactive{} && tar -zxf {}.tar.gz'.format(self.executionpath, self.interactive, dirname)
        else:
            if not isempty(restart):
                launchcommand = 'cd {} && cd {} && qsub {}.queue'.format(self.executionpath, dirname, modelname)
            else:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz  && qsub {}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)

        #Execute Queue job
        issmssh(self.name, self.login, self.port, launchcommand)

    # }}}
    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        directory = '{}/{}/'.format(self.executionpath, dirname)
        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
