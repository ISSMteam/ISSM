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
        self.name           = 'gadi.nci.org.au'
        self.login          = ''
        # Adjust modules for Gadi as needed:
        self.modules        = [
            # Example modules or spack loads:
            'openmpi/4.1.2', 
            # 'python3/3.9.2',
            # or just keep them empty if you handle externally:
        ]
        # For Gadi, you usually specify the total ncpus. But if you still want
        # to specify by numnodes & cpuspernode, that’s fine. Adjust as needed:
        self.numnodes       = 1
        self.cpuspernode    = 4
        self.port           = ''  # typical SSH port
        self.queue          = 'normal'  # or "hugemem", "express", etc.
        self.time           = 60  # total minutes of walltime, e.g. 60 => 1 hour
        self.processor      = ''  # not usually needed for Gadi
        self.srcpath        = ''
        self.extpkgpath     = ''
        self.codepath       = ''
        self.executionpath  = ''
        self.project        = ''  # Gadi uses -P <PROJECT>
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
        s = 'class Gadi object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    modules: {}\n'.format(strjoin(self.modules, ', '))
        s += '    numnodes: {}\n'.format(self.numnodes)
        s += '    cpuspernode: {}\n'.format(self.cpuspernode)
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
        return self
    # }}}

    def BuildQueueScript(
        self, dirname, modelname, solution,
        io_gather, isvalgrind, isgprof, isdakota, isoceancoupling
    ):  # {{{
        """
        Create a PBS script for Gadi. 
        Gadi typically uses #PBS lines like:
         - #PBS -P <PROJECT>
         - #PBS -q normal
         - #PBS -l ncpus=...,walltime=HH:MM:SS,mem=...
         - #PBS -l wd
         - #PBS -j oe
        """

        if isgprof:
            print('gprof not typically used on Gadi via this script, ignoring...')

        executable = 'issm.exe'
        if isdakota:
            version_str = IssmConfig('_DAKOTA_VERSION_')
            # e.g. "6.14"
            version_float = float(version_str[0:3])
            if version_float >= 6.0:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Convert self.time (minutes) to hh:mm:ss
        hours   = self.time // 60
        minutes = self.time % 60
        walltime_str = '{:02d}:{:02d}:00'.format(hours, minutes)

        # Write queue script
        fid = open(modelname + '.queue', 'w')
        fid.write('#!/bin/bash\n')
        fid.write('#PBS -P {}\n'.format(self.project))
        fid.write('#PBS -q {}\n'.format(self.queue))
        fid.write('#PBS -l ncpus={}\n'.format(self.nprocs()))
        # Optionally request memory (example: 4GB x nprocs):
        # fid.write('#PBS -l mem={}GB\n'.format(4 * self.nprocs()))
        fid.write('#PBS -l walltime={}\n'.format(walltime_str))
        fid.write('#PBS -l wd\n')  
        fid.write('#PBS -j oe\n')
        fid.write('#PBS -l storage=scratch/{0}+gdata/{0}\n'.format(self.project))
        fid.write('#PBS -m bea\n')
        fid.write('#PBS -o {}/{}/{}.outlog \n'.format(self.executionpath, dirname, modelname))
        fid.write('#PBS -e {}/{}/{}.errlog \n\n'.format(self.executionpath, dirname, modelname))   

        fid.write('\n# Load modules as needed:\n')
        fid.write('module purge\n')
        for m in self.modules:
            fid.write('module load {}\n'.format(m))

        # Optionally source environment scripts if needed:
        # fid.write('source /g/data/...someSpackOrCondaEnv...\n')

        fid.write('\n# Switch to run directory (if not using -l wd):\n')
        fid.write('cd {}/{}/\n'.format(self.executionpath, dirname))
        fid.write('\n# Now launch the job:\n')
        

        if not isvalgrind:
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
        # Compress inputs into a tarball
        compressstring = 'tar -zcf {}.tar.gz'.format(dirname)
        for f in filelist:
            compressstring += ' {}'.format(f)

        subprocess.call(compressstring, shell=True)
        print('Uploading input file and queueing script to Gadi...')
        directory = self.executionpath
        issmscpout(self.name, directory, self.login, self.port, [dirname + '.tar.gz'])
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        """
        On Gadi, typically you do: 
          qsub modelname.queue
        rather than `./modelname.queue`.
        """
        if not isempty(restart):
            # If "restart" logic is needed, adapt as necessary
            launchcommand = (
                'cd {executionpath}/{dirname} && qsub {modelname}.queue'
                .format(executionpath=self.executionpath, 
                        dirname=dirname, modelname=modelname)
            )
        else:
            # Create/clean directory, extract tar, then qsub
            launchcommand = (
                'cd {executionpath} && '
                'rm -rf ./{dirname} && mkdir {dirname} && cd {dirname} && '
                'mv ../{dirname}.tar.gz ./ && '
                'tar -zxf {dirname}.tar.gz && '
                'qsub {modelname}.queue'
                .format(executionpath=self.executionpath, 
                        dirname=dirname, modelname=modelname)
            )
            
        print('Launching solution sequence on Gadi via SSH...')
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from Gadi back to local machine
        directory = '{}/{}/'.format(self.executionpath, dirname)
        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
