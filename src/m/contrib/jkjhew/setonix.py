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
    from setonix_settings import setonix_settings
except ImportError:
    pass
from QueueRequirements import QueueRequirements

class setonix(object):
    """Setonix cluster class definition (Pawsey).
       Usage:
           cluster = setonix()  # default constructor
           cluster = setonix('numnodes',2,'cpuspernode',64,'login','username',...)
    """

    def __init__(self, *args):  # {{{
        self.name           = 'setonix.pawsey.org.au'
        self.login          = ''
        # Adjust modules as needed for your Setonix environment:
        self.modules        = [
            # e.g. 'cpeCray', 'openmpi', or any relevant modules
            # 'cpeGNU/22.08', 'openmpi/4.1.4', etc.
        ]
        # For Setonix/Slurm, you typically specify total cores via numnodes*cpuspernode,
        # or just set them directly.  Adjust as needed:
        self.numnodes       = 1
        self.cpuspernode    = 4
        self.port           = ''  
        # Slurm uses "partitions" instead of "queues". We'll keep the field named "queue"
        # to stay consistent with other cluster classes, but you can rename if you like.
        self.queue          = 'work'  
        self.time           = 60  # total minutes of walltime
        self.processor      = ''  # not always necessary
        self.srcpath        = ''
        self.extpkgpath     = ''
        self.codepath       = ''
        self.executionpath  = ''
        self.project        = ''  # Slurm account => e.g. "--account=xx00"
        self.interactive    = 0
        self.bbftp          = 0
        self.numstreams     = 8
        self.hyperthreading = 0

        # Use provided options to change fields
        options = pairoptions(*args)

        # Try user settings file if you have one
        try:
            self = setonix_settings(self)
        except NameError:
            pass

        # Override defaults with user-specified parameters
        self = options.AssignObjectFields(self)

        # Ensure login is set correctly
        if not self.login:
            raise ValueError("Login must be specified for Setonix cluster.")
    # }}}

    def __repr__(self):  # {{{
        s = 'class Setonix object\n'
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    modules: {}\n'.format(strjoin(self.modules, ', '))
        s += '    numnodes: {}\n'.format(self.numnodes)
        s += '    cpuspernode: {}\n'.format(self.cpuspernode)
        s += '    np: {}\n'.format(self.nprocs())
        s += '    port: {}\n'.format(self.port)
        s += '    queue (partition): {}\n'.format(self.queue)
        s += '    time (min): {}\n'.format(self.time)
        s += '    processor: {}\n'.format(self.processor)
        s += '    srcpath: {}\n'.format(self.srcpath)
        s += '    extpkgpath: {}\n'.format(self.extpkgpath)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    project (account): {}\n'.format(self.project)
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
        # Example: define some partitions and their max time/proc
        # (TOTALLY EXAMPLE – adapt to your actual Pawsey constraints!)
        queuedict = {
            'work': [48*60,  1024],  # up to 48h, 1024 total cores
            'debug': [2*60, 128],    # up to 2h, 128 cores
        }
        QueueRequirements(queuedict, self.queue, self.nprocs(), self.time)

        # Basic checks
        if not self.login:
            md = md.checkmessage('login is empty')
        if not self.srcpath:
            md = md.checkmessage('srcpath is empty')
        if not self.codepath:
            md = md.checkmessage('codepath is empty')
        if not self.executionpath:
            md = md.checkmessage('executionpath is empty')
        if not self.project:
            md = md.checkmessage('project (Slurm account) is empty')
        return self
    # }}}

    def BuildQueueScript(
        self, dirname, modelname, solution,
        io_gather, isvalgrind, isgprof, isdakota, isoceancoupling
    ):  # {{{
        """
        Create a Slurm script for Setonix (Pawsey).
        Typically uses lines like:
         #SBATCH --account=<PROJECT>
         #SBATCH --partition=<QUEUE>
         #SBATCH --time=HH:MM:SS
         #SBATCH --nodes=<NUMNODES>
         #SBATCH --ntasks=<TOTALTASKS>
         #SBATCH --ntasks-per-node=<TASKSPERNODE>
         ...
        """

        if isgprof:
            print('gprof not typically invoked this way, ignoring...')

        # Decide which binary to run
        executable = 'issm.exe'
        if isdakota:
            version_str = IssmConfig('_DAKOTA_VERSION_')
            version_float = float(version_str[0:3])  # e.g. "6.14" -> 6.1
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
        fid.write('#\n')
        fid.write('# Slurm job script for Setonix (Pawsey)\n\n')

        fid.write('#SBATCH --account={}\n'.format(self.project))
        fid.write('#SBATCH --partition={}\n'.format(self.queue))
        fid.write('#SBATCH --time={}\n'.format(walltime_str))
        fid.write('#SBATCH --nodes={}\n'.format(self.numnodes))
        fid.write('#SBATCH --ntasks={}\n'.format(self.nprocs()))
        fid.write('#SBATCH --ntasks-per-node={}\n'.format(self.cpuspernode))
        # If you need memory or specific job name, add:
        # fid.write('#SBATCH --mem=0\n')  # Just an example, "0" can be "no limit" on some systems
        fid.write('#SBATCH --job-name={}\n'.format(modelname))
        fid.write('#SBATCH --output={}/{}/{}.outlog\n'.format(self.executionpath, dirname, modelname))
        fid.write('#SBATCH --error={}/{}/{}.errlog\n'.format(self.executionpath, dirname, modelname))
        fid.write('\n')

        fid.write('echo "Starting job on `date`"\n\n')

        fid.write('# Load modules as needed:\n')
        fid.write('module purge\n')
        for m in self.modules:
            fid.write('module load {}\n'.format(m))

        # Typically on Pawsey, you might need to load the Cray compiler env or MPI, e.g.:
        # fid.write('module load cpeGNU/22.08\n')
        # fid.write('module load cray-mpich\n')
        # Adjust to your environment or if you handle modules outside this script, remove.

        fid.write('\n# Move to run directory:\n')
        fid.write('cd {}/{}/\n'.format(self.executionpath, dirname))

        fid.write('\n# Launch the job:\n')
        # Use srun (Slurm's recommended launching method)
        if not isvalgrind:
            fid.write('srun {}/{} {} {}/{} {}\n'.format(
                self.codepath, executable,
                solution, self.executionpath, dirname, modelname))
        else:
            # Example valgrind usage
            supstring = ''
            # if you have self.valgrindsup, loop here:
            # for supfile in self.valgrindsup:
            #     supstring += ' --suppressions=' + supfile
            fid.write(
                'srun valgrind --leak-check=full{} {}/{} {} {}/{} {}\n'.format(
                    supstring, self.codepath, executable,
                    solution, self.executionpath, dirname, modelname))

        # If your run generates multiple *.outbin.# partial files and
        # you need to combine them, add something like:
        if not io_gather:
            fid.write('\ncat {}.outbin.* > {}.outbin\n'.format(modelname, modelname))

        fid.write('\necho "Job finished on `date`"\n')
        fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress local files into a tarball
        compressstring = 'tar -zcf {}.tar.gz'.format(dirname)
        for f in filelist:
            compressstring += ' {}'.format(f)
        subprocess.call(compressstring, shell=True)

        print('Uploading input files and queue script to Setonix...')
        directory = self.executionpath
        issmscpout(self.name, directory, self.login, self.port, [dirname + '.tar.gz'])
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        """
        On Setonix, typically we do:
          sbatch modelname.queue
        rather than something like "./modelname.queue".
        """

        if not isempty(restart):
            # If "restart" logic is needed, adapt as necessary
            launchcommand = (
                'cd {executionpath}/{dirname} && sbatch {modelname}.queue'
                .format(executionpath=self.executionpath,
                        dirname=dirname, modelname=modelname)
            )
        else:
            # Create/clean directory, extract tar, then sbatch
            launchcommand = (
                'cd {executionpath} && '
                'rm -rf ./{dirname} && mkdir {dirname} && cd {dirname} && '
                'mv ../{dirname}.tar.gz ./ && '
                'tar -zxf {dirname}.tar.gz && '
                'sbatch {modelname}.queue'
                .format(executionpath=self.executionpath,
                        dirname=dirname, modelname=modelname)
            )

        print('Launching solution sequence on Setonix via SSH...')
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files back from Setonix to local machine
        directory = '{}/{}/'.format(self.executionpath, dirname)
        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
