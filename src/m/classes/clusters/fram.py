import subprocess

import numpy as np

from fielddisplay import fielddisplay
try:
    from fram_settings import fram_settings
except ImportError:
    print('You need fram_settings.py to proceed, check presence and sys.path')
from helpers import *
from pairoptions import pairoptions
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from QueueRequirements import QueueRequirements


class fram(object):
    """FRAM cluster class definition

    This is a SLURM queue
    The priorities are based on a point system, reservation when reaching 20000 and earning 1 point per min.
    - Devel queue starts at 19990
    - Normal starts at 19940
    - Normal unpri atarts at 19400

    Jobs can be:
    - normal (4 to 30 nodes, more if asked, 48h max walltime, 60Gb per nodes)
    - bigmem for big memory nodes (8 512Gb nodes and 2 6Tb nodes, shared nodes, 14days max walltime

    Usage:
        cluster = fram()
    """

    def __init__(self, *args):  # {{{
        self.name = 'fram'
        self.login = ''
        self.numnodes = 2
        self.cpuspernode = 20
        self.mem = 1.6
        self.queue = 'normal'
        self.time = 2 * 60
        self.codepath = ''
        self.executionpath = ''
        self.interactive = 0
        self.port = []
        self.accountname = ''
        self.profiling = 0
        # Use provided options to change fields
        options = pairoptions(*args)

        # Initialize cluster using user settings if provided
        self = fram_settings(self)

        # OK get other fields
        self = options.AssignObjectFields(self)
        self.np = self.numnodes * self.cpuspernode
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = "class fram object:"
        s = "%s\n%s" % (s, fielddisplay(self, 'name', 'name of the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'login', 'login'))
        s = "%s\n%s" % (s, fielddisplay(self, 'numnodes', 'number of nodes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'cpuspernode', 'number of nodes per CPUs'))
        s = "%s\n%s" % (s, fielddisplay(self, 'mem', 'memory per CPU'))
        s = "%s\n%s" % (s, fielddisplay(self, 'queue', 'name of the queue (normal (D), short, singlenode, multinode, devel)'))
        s = "%s\n%s" % (s, fielddisplay(self, 'time', 'walltime requested in minutes'))
        s = "%s\n%s" % (s, fielddisplay(self, 'codepath', 'code path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'executionpath', 'execution path on the cluster'))
        s = "%s\n%s" % (s, fielddisplay(self, 'interactive', ''))
        s = "%s\n%s" % (s, fielddisplay(self, 'accountname', 'your cluster account'))
        s = "%s\n%s" % (s, fielddisplay(self, 'profiling', 'enable profiling if 1 default is 0'))
        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Queue dictionary gives queue name as key and max walltime and CPUs as var
        queuedict = {'normal': [2 * 24 * 60, 2048],
                     'devel': [4 * 60, 2048]}
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
        fid = open(modelname + '.queue', 'w')

        fid.write('#!/bin/bash -l\n')
        fid.write('#SBATCH --job-name=%s \n' % shortname)
        fid.write('#SBATCH --partition %s \n' % self.queue)
        fid.write('#SBATCH --nodes=%i' % self.numnodes)
        fid.write('#SBATCH --ntasks-per-nodes==%i \n' % self.cpuspernode)
        fid.write('#SBATCH --time=%s\n' % self.time)  #walltime is minutes
        fid.write('#SBATCH --mem-per-cpu=%iGB\n' % self.mem)  # mem is in GB
        if (np.mod(self.np, 16) + np.mod(self.np, 20)) == 0:
            fid.write('#SBATCH --ntask=%i\n' % self.np)
        fid.write('#SBATCH --account=%s\n' % self.accountname)
        fid.write('#SBATCH --output %s/%s /%s.outlog \n' % (self.executionpath, dirname, modelname))
        fid.write('#SBATCH --error %s/%s /%s.errlog \n\n' % (self.executionpath, dirname, modelname))

        fid.write('export ISSM_DIR="%s/../ "\n' % self.codepath)
        fid.write('module restore system\n')
        fid.write('module load load Automake/1.15.1-GCCcore-6.3.0\n')
        fid.write('module load libtool/2.4.6-GCCcore-6.3.0\n')
        fid.write('module load CMake/3.9.1\n')
        fid.write('module load PETSc/3.8.0-intel-2017a-Python-2.7.13\n')
        fid.write('module load ParMETIS/4.0.3-intel-2017a\n')
        fid.write('cd %s/%s/ \n\n' % (self.executionpath, dirname))
        if self.profiling:
            fid.write('module load perf-report\n')
            fid.write('perf-report mpirun -np %i %s/%s %s %s/%s %s\n' % (self.np, self.codepath, executable, str(solution), self.executionpath, dirname, modelname))
        else:
            fid.write('mpirun -np %i %s/%s %s %s/%s %s\n' % (self.np, self.codepath, executable, str(solution), self.executionpath, dirname, modelname))
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
        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
