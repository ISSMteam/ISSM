from subprocess import call

import numpy as np

from IssmConfig import IssmConfig
from issmdir import *
from pairoptions import pairoptions


class local(object):  # {{{
    """LOCAL class definition

    Usage:
        cluster = local('name', 'astrid', 'np', 3)
        cluster = local('name', oshostname(), 'np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name           = ''
        self.np             = 1
        self.codepath       = IssmConfig('ISSM_PREFIX') + '/bin'
        self.etcpath        = issmdir() + '/etc'
        self.executionpath  = issmdir() + '/execution'
        self.verbose        = 1
        self.shell          = '/bin/sh'

        # Use provided options to change fields
        options = pairoptions(*args)

        # Get name
        cluster.name = options.getfieldvalue('name', oshostname())

        # OK get other fields
        self = options.AssignObjectFields(self)

    def __repr__(cluster):  # {{{
        # Display the object
        s = 'class {} = \n'.format(type(cluster).__name__)
        s += '    name: {}\n'.format(cluster.name)
        s += '    np: {}\n'.format(cluster.np)
        s += '    codepath: {}\n'.format(cluster.codepath)
        s += '    executionpath: {}\n'.format(cluster.executionpath)
        s += '    etcpath: {}\n'.format(cluster.etcpath)
        s += '    verbose: {}\n'.format(cluster.verbose)
        s += '    shell: {}\n'.format(cluster.shell)

        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if cluster.np < 1:
            md.checkmessage('number of processors should be at least 1')

        if np.isnan(cluster.np):
            md.checkmessage('number of processors should not be NaN')

        return md
    # }}}

    def BuildQueueScript(cluster, dirname, modelname, solution, io_gather, isvalgrind, isgporf, isdakota, isoceancoupling):  # {{{
        # Which executable are we calling?
        executable = 'issm.exe' # Default

        if isdakota:
            executable = 'issm_dakota.exe'

        fid = open(modelname + '.queue', 'w')
        fid.write('#!{}\n'.format(cluster.shell))
        fid.write('mpiexec -np {} {}/{} {} {} {}\n',cluster.np,cluster.codepath,executable,solution,'./',modelname)
        fid.close()

    def UploadQueueJob(cluster, modelname, dirname, filelist):  # {{{
        # Do nothing really
        pass
    # }}}

    def LaunchQueueJob(cluster, modelname, dirname, filelist, restart, batch):  # {{{
        subprocess.call('source ' + modelname + '.queue')
    # }}}

    def Download(cluster, dirname, filelist):  # {{{
        pass
    # }}}
