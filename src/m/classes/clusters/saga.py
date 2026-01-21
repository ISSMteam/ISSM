import datetime
import subprocess

from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from pairoptions import pairoptions
from QueueRequirements import QueueRequirements
try:
    from saga_settings import saga_settings
except ImportError:
    print('You need saga_settings.py to proceed, check presence and sys.path')


class saga(object):
    """SAGA cluster class definition

    This is a SLURM queue

    Usage:
        cluster = saga()
    """

    def __init__(self, *args):  # {{{
        self.name = 'saga'

        self.numnodes = 1
        self.cpuspernode = 20
        self.mem = 2
        self.queue = 'normal'
        self.time = 2 * 60

        self.port = []
        #set by setting file
        self.login = ''
        self.codepath = ''
        self.executionpath = ''
        self.accountname = ''

        self.interactive = 0
        self.profiling = 0

        self.valgrind = '/cluster/software/Valgrind/3.16.1-gompi-2019b/bin/valgrind'

        # Use provided options to change fields
        options = pairoptions(*args)
        # Initialize cluster using user settings if provided
        self = saga_settings(self)
        # OK get other fields
        self = options.AssignObjectFields(self)
        self.np = self.numnodes * self.cpuspernode
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = "class vilje object:"
        s = "%s\n%s" % (s, fielddisplay(self, 'name', 'name of the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'login', 'login'))
        s = "%s\n%s" % (s, fielddisplay(self, 'numnodes', 'number of nodes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'cpuspernode', 'number of CPUs per nodes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'mem', 'memory per CPU'))
        s = "%s\n%s" % (s, fielddisplay(self, 'queue', 'name of the queue (normal (D), bigmem, devel)'))
        s = "%s\n%s" % (s, fielddisplay(self, 'time', 'walltime requested in minutes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'codepath', 'code path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'executionpath', 'execution path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'interactive', ''))
        s = "%s\n%s" % (s, fielddisplay(self, 'accountname', 'your cluster account'))
        s = "%s\n%s" % (s, fielddisplay(self, 'profiling', 'enable profiling if 1 default is 0'))
        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Queue dictionarry  gives queue name as key and max walltime, CPUs, and memory (GB) as var
        queuedict = {'normal': [7 * 24 * 60, 256, 8],
                     'bigmem': [14 * 24 * 60, 256, 10],
                     'devel': [2 * 60, 256, 8]}
        QueueRequirements(queuedict, self.queue, self.np, self.time)

        # Miscellaneous
        if not self.login:
            md = md.checkmessage('login empty')
        if not self.codepath:
            md = md.checkmessage('codepath empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath empty')
        if self.interactive == 1:
            md = md.checkmessage('interactive mode not implemented')
        return self
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        executable = 'issm.exe'
        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')[0:2]
            version = float(version)
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'
        # Write queuing script
        shortname = modelname[0:min(12, len(modelname))]
        timeobj = datetime.timedelta(minutes=self.time)
        m, s = divmod(timeobj.total_seconds(), 60)
        h, m = divmod(m, 60)
        d, h = divmod(h, 24)
        timestring = "%02d-%02d:%02d:%02d" % (d, h, m, s)
        print('timestring')
        fid = open(modelname + '.queue', 'w')
        fid.write('#!/bin/bash -l\n')
        fid.write('#SBATCH --job-name=%s \n' % shortname)
        if self.queue in ['devel']:
            fid.write('#SBATCH --partition=normal \n')
            fid.write('#SBATCH --qos=%s \n' % self.queue)
        else:
            fid.write('#SBATCH --partition=%s \n' % self.queue)

        fid.write('#SBATCH --nodes=%i \n' % self.numnodes)
        fid.write('#SBATCH --ntasks=%i \n' % self.cpuspernode)
        fid.write('#SBATCH --time={}\n'.format(timestring))  #walltime is minutes
        fid.write('#SBATCH --mem-per-cpu={}M\n'.format(int(1000 * self.mem)))  # mem is in MB

        fid.write('#SBATCH --account=%s\n' % self.accountname)
        fid.write('#SBATCH --output %s/%s/%s.outlog \n' % (self.executionpath, dirname, modelname))
        fid.write('#SBATCH --error %s/%s/%s.errlog \n\n' % (self.executionpath, dirname, modelname))

        fid.write('export ISSM_DIR="%s/../"\n' % self.codepath)
        fid.write('module purge\n')
        fid.write('module load CMake/3.15.3-GCCcore-8.3.0\n')
        fid.write('module load PETSc/3.12.4-foss-2019b\n')
        fid.write('module load ParMETIS/4.0.3-gompi-2019b\n')
        if isvalgrind:
            fid.write('module --ignore-cache load Valgrind/3.16.1-gompi-2019b \n')

        fid.write('cd %s/%s/ \n\n' % (self.executionpath, dirname))
        if isvalgrind:
            # profiling
            #fid.write('srun {} --tool=callgrind {}/{} {} {}/{} {} 2>{}.errlog>{}.outlog \n'.format(self.valgrind, self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
            # leak check
            fid.write('mpirun --bind-to none {} --leak-check=full {}/{} {} {}/{} {} 2>{}.errlog>{}.outlog '.format(self.valgrind, self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
        else:
            fid.write('time mpirun --bind-to none {}/{} {} {}/{} {}\n'.format(self.codepath, executable, solution, self.executionpath, dirname, modelname))
        fid.close()

    # }}}
    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one zip
        compressstring = 'tar -zcf %s.tar.gz ' % dirname
        for file in filelist:
            compressstring += ' {}'.format(file)
        subprocess.call(compressstring, shell=True)

        #upload input files
        issmscpout(self.name, self.executionpath, self.login, self.port, [dirname + '.tar.gz'])
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        #Execute Queue job
        if not isempty(restart):
            launchcommand = 'cd %s && cd %s && sbatch %s.queue' % (self.executionpath, dirname, modelname)
        else:
            launchcommand = 'cd %s && rm -rf ./%s && mkdir %s && cd %s && mv ../%s.tar.gz ./ && tar -zxf %s.tar.gz  && sbatch %s.queue' % (self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        directory = '%s/%s/' % (self.executionpath, dirname)
        try:
            issmscpin(self.name, self.login, self.port, directory, filelist)
        except OSError:
            print("File does not exsit, skiping")
            # }}}
