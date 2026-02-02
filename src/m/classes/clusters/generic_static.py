import math
import os
import subprocess
import numpy as np
from IssmConfig import IssmConfig
from issmdir import issmdir
from issmssh import issmssh
from issmscpin import issmscpin
from issmscpout import issmscpout
from MatlabFuncs import *
from pairoptions import pairoptions


class generic_static(object):
    """GENERIC cluster class definition

    Usage:
        cluster = generic_static('name', 'astrid', 'np', 3)
    """

    def __init__(self, *args):  # {{{
        codepath = subprocess.check_output(["which", "issm.exe"]).rstrip('\r\n')
        codepath = codepath.replace('/issm.exe', '')

        self.name = ''
        self.np = 1
        self.codepath = codepath
        self.executionpath = '.'
        self.interactive = 1
        self.shell = '/bin/sh'

        #use provided options to change fields
        options = pairoptions(*args)

        #get name
        self.name = oshostname()

        #initialize cluster using user settings if provided
        if os.path.exists(self.name + '_settings.py'):
            exec(compile(open(self.name + '_settings.py').read(), self.name + '_settings.py', 'exec'), globals())

        #OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        #  display the object
        s = "class '%s' object '%s' = \n" % (type(self), 'self')
        s += "    name: %s\n" % self.name
        s += "    np: %i\n" % self.np
        s += "    codepath: %s\n" % self.codepath
        s += "    shell: %s\n" % self.shell
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.np
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if self.np < 1:
            md.checkmessage('number of processors should be at least 1')
        if np.isnan(self.np):
            md.checkmessage('number of processors should not be NaN!')

        return md
    # }}}
    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        # Which executable are we calling?
        executable = 'issm.exe' # default

        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')
            version = float(version[0])
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Check that executable exists at the right path
        if not os.path.isfile(self.codepath + '/' + executable):
            raise RuntimeError('File ' + self.codepath + '/' + executable + ' does not exist')

        # Process codepath and prepend empty spaces with \ to avoid errors in queuing script
        codepath = self.codepath.replace(' ', r'\ ')

        # Write queuing script
        fid = open(modelname + '.queue', 'w')
        fid.write('#!{}'.format(self.shell) + '\n')
        fid.write('{}/mpiexec -np {} {}/{} {} {} {}'.format(codepath, self.np, codepath, executable, solution, './', modelname))
        fid.close()

        # Set permissions on queue script so that it can be run
        subprocess.call(['chmod', '0755', modelname + '.queue'])

        # Create an errlog and outlog file
        fid = open(modelname + '.errlog', 'w')
        fid.close()
        fid = open(modelname + '.outlog', 'w')
        fid.close()
    # }}}

    def BuildKrigingQueueScript(self, modelname, solution, io_gather, isvalgrind, isgprof):  # {{{
        # Which executable are we calling?
        executable = 'kriging.exe' # default

        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')
            version = float(version[0])
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Check that executable exists at the right path
        if not os.path.isfile(self.codepath + '/' + executable):
            raise RuntimeError('File ' + self.codepath + '/' + executable + ' does not exist')

        # Process codepath and prepend empty spaces with \ to avoid errors in queuing script
        codepath = self.codepath.replace(' ', r'\ ')

        # Write queuing script
        fid = open(modelname + '.queue', 'w')
        fid.write('#!{}'.format(self.shell) + '\n')
        fid.write('{}/mpiexec -np {} {}/{} {} {} {}'.format(codepath, self.np, codepath, executable, solution, './', modelname) + '\n')
        fid.close()

        # Set permissions on queue script so that it can be run
        subprocess.call(['chmod', '0755', modelname + '.queue'])

        # Create an errlog and outlog file
        fid = open(modelname + '.errlog', 'w')
        fid.close()
        fid = open(modelname + '.outlog', 'w')
        fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Do nothing
        return
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        if not ispc():
            # Figure out which file extension to use
            if self.shell.find('csh') == -1:
                shellext='sh'
            else:
                shellext='csh'

            launchcommand = './' + modelname + '.queue'
            subprocess.call([launchcommand])
        else:
            launchcommand = './' + modelname + '.bat'
            subprocess.call([launchcommand])
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Do nothing
        return
    # }}}
