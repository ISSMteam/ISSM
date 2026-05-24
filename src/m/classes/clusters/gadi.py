import os
import subprocess

import cluster_defaults
from fielddisplay import fielddisplay
from helpers import *
from IssmConfig import IssmConfig
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from issmdir import issmdir
from MatlabFuncs import *
from pairoptions import pairoptions
try:
    from gadi_settings import gadi_settings
except ImportError:
    pass
from QueueRequirements import QueueRequirements

class gadi(object):
    """Gadi cluster class definition.
       Usage:
           cluster = gadi()  # default constructor
           cluster = gadi('numnodes',2,'cpuspernode',24,'login','username') ...
    """

    def __init__(self, *args):  # {{{
        self.name           = oshostname()
        self.login          = ''
        self.moduleload     = []
        self.moduleuse      = []
        self.numnodes       = 1
        self.cpuspernode    = 4
        self.memory         = 40  # e.g. '40GB'
        self.port           = 0   # typical SSH port
        self.queue          = 'normal'  # or "hugemem", "express", etc.
        self.time           = 60  # total minutes of walltime, e.g. 60 => 1 hour
        self.processor      = ''  # not usually needed for Gadi
        self.srcpath        = ''
        self.extpkgpath     = ''
        self.codepath       = ''
        self.executionpath  = ''
        self.project        = ''
        self.storage        = ''
        self.interactive    = 0
        self.bbftp          = 0
        self.numstreams     = 8
        self.hyperthreading = 0

        # Use provided options to change fields
        options = pairoptions(*args)

        # If you had a user settings file, try it:
        try:
            self = gadi_settings(self)
        except NameError:
            pass

        # Override defaults with user-specified parameters
        self = options.AssignObjectFields(self)

        # Ensure login is set correctly
        if not self.login:
            raise ValueError("Login must be specified for Gadi cluster.")
        
    # }}}

    def __repr__(self):  # {{{
        s = 'class gadi object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    moduleuse: {}\n'.format(', '.join(self.moduleuse) if getattr(self, 'moduleuse', None) else '')
        s += '    moduleload: {}\n'.format(', '.join(self.moduleload) if getattr(self, 'moduleload', None) else '')
        s += '    numnodes: {}\n'.format(self.numnodes)
        s += '    cpuspernode: {}\n'.format(self.cpuspernode)
        s += '    memory: {}GB\n'.format(self.memory)
        s += '    np: {}\n'.format(self.nprocs())
        s += '    port: {}\n'.format(self.port)
        s += '    queue: {}\n'.format(self.queue)
        s += '    time (min): {}\n'.format(self.time)
        s += '    processor: {}\n'.format(self.processor)
        s += '    srcpath: {}\n'.format(self.srcpath)
        s += '    extpkgpath: {}\n'.format(self.extpkgpath)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    project: {}\n'.format(self.project)
        s += '    storage: {}\n'.format(self.storage)
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
        queuedict = {
            'normal': [48*60,  3072],  # e.g. up to 48h, 3072 cores
            'express': [2*60,  960],   # e.g. up to 2h, 960 cores
            'hugemem': [48*60, 3072],
        }
        QueueRequirements(queuedict, self.queue, self.nprocs(), self.time)

        # Some minimal checks
        if not self.login:
            md = md.checkmessage('login is empty')
        if not self.srcpath:
            md = md.checkmessage('srcpath is empty')
        if not self.codepath:
            md = md.checkmessage('codepath is empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath is empty')
        if not self.project:
            md = md.checkmessage('project (PBS -P) is empty')
        if not self.storage:
            md = md.checkmessage('storage is empty')
        if len(self.moduleload) != len(self.moduleuse):
            raise ValueError("moduleload and moduleuse must have the same length")
        return self
    # }}}

    def BuildQueueScript(self, md, filename, executable):  # {{{
        """
        Create a PBS script for Gadi. 
        Gadi typically uses #PBS lines like:
         - #PBS -P <PROJECT>
         - #PBS -q normal
         - #PBS -l ncpus=...,walltime=HH:MM:SS,mem=...
         - #PBS -l wd
         - #PBS -j oe
        """

        # Get variables from md
        dirname         = md.private.runtimename
        modelname       = md.miscellaneous.name
        solution        = md.private.solution
        io_gather       = md.settings.io_gather

        if md.debug.gprof:
            print('gprof not typically used on Gadi via this script, ignoring...')

        # Convert self.time (minutes) to hh:mm:ss
        hours   = self.time // 60
        minutes = self.time % 60
        walltime_str = '{:02d}:{:02d}:00'.format(hours, minutes)

        # Write queue script
        fid = open(filename, 'w')
        fid.write('#!/bin/bash\n')
        fid.write('#PBS -P {}\n'.format(self.project))
        fid.write('#PBS -q {}\n'.format(self.queue))
        fid.write('#PBS -l ncpus={}\n'.format(self.nprocs()))
        fid.write('#PBS -l mem={}GB\n'.format(self.memory))
        fid.write('#PBS -l walltime={}\n'.format(walltime_str))
        fid.write('#PBS -l wd\n')  
        fid.write('#PBS -j oe\n')
        fid.write('#PBS -l storage={}\n'.format(self.storage))
        fid.write('#PBS -m bea\n')
        fid.write('#PBS -o {}/{}/{}.outlog \n'.format(self.executionpath, dirname, modelname))
        fid.write('#PBS -e {}/{}/{}.errlog \n\n'.format(self.executionpath, dirname, modelname))   

        fid.write('\n# Load modules as needed:\n')
        fid.write('module purge\n')
        # Print alternating lines
        for x, y in zip(self.moduleuse, self.moduleload):
            fid.write(f"module use {x}\n")
            fid.write(f"module load {y}\n")

        # Optionally source environment scripts if needed:
        # fid.write('source /g/data/...someSpackOrCondaEnv...\n')

        fid.write('\n# Switch to run directory (if not using -l wd):\n')
        fid.write('cd {}/{}/\n'.format(self.executionpath, dirname))
        fid.write('\n# Now launch the job:\n')
        

        if not md.debug.valgrind:
            fid.write('mpiexec -n {} {}/{} {} {}/{} {}\n'.format(
                self.nprocs(), self.codepath, executable,
                solution, self.executionpath, dirname, modelname))
        else:
            # Example valgrind usage
            supstring = ''
            # If you have self.valgrindsup or so, otherwise remove
            # for supfile in self.valgrindsup:
            #     supstring += ' --suppressions=' + supfile
            fid.write(
                'mpiexec -n {} valgrind --leak-check=full{} {}/{} {} {}/{} {}\n'.format(
                    self.nprocs(), supstring, self.codepath, executable,
                    solution, self.executionpath, dirname, modelname))

        if not io_gather:
            # Merge partial outbin files if your code writes multiple .outbin.#
            fid.write('\ncat {}.outbin.* > {}.outbin\n'.format(modelname, modelname))

        fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        cluster_defaults.UploadQueueJob(self, modelname, dirname, filelist)
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        cluster_defaults.LaunchQueueJobSbatch(self, modelname, dirname, filelist, restart, batch, 3)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        cluster_defaults.Download(self, dirname, filelist)
    # }}}
