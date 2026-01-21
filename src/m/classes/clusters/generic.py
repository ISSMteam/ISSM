from subprocess import call

from fielddisplay import fielddisplay
from helpers import *
try:
    from generic_settings import generic_settings
except ImportError:
    #print('Warning: generic.py: generic_settings.py not found, default will be used')
    pass
from MatlabFuncs import *
from IssmConfig import IssmConfig
from issmdir import issmdir
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh
from pairoptions import pairoptions


class generic(object):
    """generic cluster class definition

    Usage:
        cluster = generic('name', 'astrid', 'np', 3)
        cluster = generic('name', gethostname(), 'np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = ''
        self.login = ''
        self.np = 1
        self.port = 0
        self.interactive = 1
        self.codepath = IssmConfig('ISSM_PREFIX')[0] + '/bin'
        self.executionpath = issmdir() + '/execution'
        self.valgrind = issmdir() + '/externalpackages/valgrind/install/bin/valgrind'
        self.valgrindlib = issmdir() + '/externalpackages/valgrind/install/lib/libmpidebug.so'
        self.valgrindsup = [issmdir() + '/externalpackages/valgrind/issm.supp']  # add any .supp in list form as needed
        self.verbose = 1
        self.shell = '/bin/sh'

        # Use provided options to change fields
        options = pairoptions(*args)

        # Get name
        self.name = oshostname()

        # Initialize cluster using user settings if provided
        try:
            self = generic_settings(self)
        except NameError:
            # print('generic_settings.py not found, using default settings')
            pass

        # OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class \'{}\' object \'{}\' = \n'.format(type(self), 'self')
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    np: {}\n'.format(self.np)
        s += '    port: {}\n'.format(self.port)
        s += '    interactive: {}\n'.format(self.interactive)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    valgrind: {}\n'.format(self.valgrind)
        s += '    valgrindlib: {}\n'.format(self.valgrindlib)
        s += '    valgrindsup: {}\n'.format(self.valgrindsup)
        s += '    verbose: {}\n'.format(self.verbose)
        s += '    shell: {}\n'.format(self.shell)
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
        executable = 'issm.exe'  # default

        if isdakota:
            version = IssmConfig('_DAKOTA_VERSION_')
            version = float(version[0])
            if version >= 6:
                executable = 'issm_dakota.exe'
        if isoceancoupling:
            executable = 'issm_ocean.exe'

        # Write queuing script
        if not ispc():
            fid = open(modelname + '.queue', 'w')
            fid.write('#!/bin/sh\n')
            if not isvalgrind:
                if self.interactive:
                    if IssmConfig('_HAVE_MPI_')[0]:
                        fid.write('mpiexec -np {} {}/{} {} {}/{} {}'.format(self.np, self.codepath, executable, solution, self.executionpath, dirname, modelname))
                    else:
                        fid.write('{}/{} {} {}/{} {}'.format(self.codepath, executable, solution, self.executionpath, dirname, modelname))
                else:
                    if IssmConfig('_HAVE_MPI_')[0]:
                        fid.write('mpiexec -np {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog'.
                                  format(self.np, self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
                    else:
                        fid.write('{}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.
                                  format(self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
            elif isgprof:
                fid.write('\n gprof {}/{} gmon.out > {}.performance'.format(self.codepath, executable, modelname))
            else:
                #Add --gen -suppressions = all to get suppression lines
                #fid.write('LD_PRELOAD={} \\\n'.format(self.valgrindlib)) it could be deleted
                supstring = ''
                for supfile in self.valgrindsup:
                    supstring += ' --suppressions=' + supfile

                if IssmConfig('_HAVE_MPI_')[0]:
                    fid.write('mpiexec -np {} {} --leak-check=full {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.
                              format(self.np, self.valgrind, supstring, self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
                else:
                    fid.write('{} --leak-check=full {} {}/{} {} {}/{} {} 2> {}.errlog > {}.outlog '.
                              format(self.valgrind, supstring, self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))

            if not io_gather:  #concatenate the output files:
                fid.write('\ncat {}.outbin .*>{}.outbin'.format(modelname, modelname))
            fid.close()

        else:  # Windows
            fid = open(modelname + '.bat', 'w')
            fid.write('@echo off\n')
            if self.interactive:
                fid.write('"{}/{}" {} "{}/{}" {} '.format(self.codepath, executable, solution, self.executionpath, dirname, modelname))
            else:
                fid.write('"{}/{}" {} "{}/{}" {} 2>{}.errlog>{}.outlog'.
                          format(self.codepath, executable, solution, self.executionpath, dirname, modelname, modelname, modelname))
            fid.close()

            #in interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(modelname + '.errlog', 'w')
            fid.close()
            fid = open(modelname + '.outlog', 'w')
            fid.close()
    # }}}

    def BuildKrigingQueueScript(self, modelname, solution, io_gather, isvalgrind, isgprof):  # {{{
        #write queuing script
        if not ispc():
            fid = open(modelname + '.queue', 'w')
            fid.write('#!/bin/sh\n')
            if not isvalgrind:
                if self.interactive:
                    fid.write('mpiexec -np {} {}/kriging.exe {}/{} {} '.format(self.np, self.codepath, self.executionpath, modelname, modelname))
                else:
                    fid.write('mpiexec -np {} {}/kriging.exe {}/{} {} 2>{}.errlog>{}.outlog '.format(self.np, self.codepath, self.executionpath, modelname, modelname, modelname, modelname))
            elif isgprof:
                fid.write('\n gprof {}/kriging.exe gmon.out>{}.performance'.format(self.codepath, modelname))
            else:
                #Add --    gen - suppressions = all to get suppression lines
                #fid.write('LD_PRELOAD={} \\\n'.format(self.valgrindlib))
                fid.write('mpiexec -np {} {} --leak -check=full --suppressions={} {}/kriging.exe {}/{} {} 2 > {}.errlog > {}.outlog ' .format(self.np, self.valgrind, self.valgrindsup, self.codepath, self.executionpath, modelname, modelname, modelname, modelname))
            if not io_gather:    #concatenate the output files:
                fid.write('\ncat {}.outbin. *>{}.outbin'.format(modelname, modelname))
            fid.close()

        else:    # Windows
            fid = open(modelname + '.bat', 'w')
            fid.write('@echo off\n')
            if self.interactive:
                fid.write('"{}/issm.exe" {} "{}/{}" {} '.format(self.codepath, solution, self.executionpath, modelname, modelname))
            else:
                fid.write('"{}/issm.exe" {} "{}/{}" {} 2>{}.errlog>{}.outlog'.format
                          (self.codepath, solution, self.executionpath, modelname, modelname, modelname, modelname))
            fid.close()

        # In interactive mode, create a run file, and errlog and outlog file
        if self.interactive:
            fid = open(modelname + '.errlog', 'w')
            fid.close()
            fid = open(modelname + '.outlog', 'w')
            fid.close()
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one zip
        compressstring = 'tar -zcf {}.tar.gz '.format(dirname)
        for file in filelist:
            compressstring += ' {}'.format(file)
        if self.interactive:
            compressstring += ' {}.errlog {}.outlog '.format(modelname, modelname)
        call(compressstring, shell=True)

        issmscpout(self.name, self.executionpath, self.login, self.port, [dirname + '.tar.gz'])

    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        if not isempty(restart):
            launchcommand = 'cd {} && cd {} chmod 755 {}.queue && ./{}.queue'.format(self.executionpath, dirname, modelname, modelname)
        else:
            if batch:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./&& tar -zxf {}.tar.gz'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname)
            else:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./&& tar -zxf {}.tar.gz  && chmod 755 {}.queue && ./{}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname, modelname)
        issmssh(self.name, self.login, self.port, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        if ispc():
            # Do nothing
            return
        # Copy files from cluster to current directory
        directory = '{}/{}/'.format(self.executionpath, dirname)
        issmscpin(self.name, self.login, self.port, directory, filelist)
    # }}}
