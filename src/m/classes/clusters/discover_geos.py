# class for running interactively on discover and writing
# binary input files for GEOS coupling
# note: need to re-implement batch script writing (i.e. .queue files) for stand-alone simulations
# that are submitted to the queue (see discover.py class for reference)
import subprocess

from helpers import *
from IssmConfig import IssmConfig
from MatlabFuncs import *
from pairoptions import pairoptions
from QueueRequirements import QueueRequirements
from generic import generic
from export_netCDF import export_netCDF
import contextlib
import io

class discover_geos(object):
    """DISCOVER cluster class definition for coupling with GEOS GCM

    Usage:
        cluster = discover()
        cluster = discover('np', 3)
        cluster = discover('np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = oshostname()
        self.login = ''
        self.modules = ['']
        self.numnodes = 1
        self.cpuspernode = 8
        self.port = 0
        self.queue = 'debug'
        self.time = 60 * 60
        self.processor = 'mil'
        self.srcpath = ''
        self.codepath = ''
        self.executionpath = ''
        self.grouplist = ''
        self.interactive = 0
        self.numstreams = 8
        self.hyperthreading = 0
        self.verbose = 1
        self.email = ''
        self.runmodel = 1       
 
        # Use provided options to change fields
        options = pairoptions(*args)

    
        # OK get other fields
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class discover_geos object\n'
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
        s += '    numstreams: {}\n'.format(self.numstreams)
        s += '    hyperthreading: {}\n'.format(self.hyperthreading)
        return s
    # }}}

    def nprocs(self):  # {{{
        return self.numnodes * self.cpuspernode
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
	# FIXME
	# version in discover.py is not up-to-date with current NCCS SCUs
        if self.nprocs() < 1:
            md.checkmessage('number of processors should be at least 1')
        if np.isnan(self.nprocs()):
            md.checkmessage('number of processors should not be NaN!')
        return md
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        # In interactive mode, create a run file, and errlog and outlog file
        executable = 'issm.exe'
        fid = open(modelname + '.run', 'w')
        fid.write('#!/bin/sh\n')
        if not isvalgrind:
            fid.write('mpiexec -np {} {}/{} {} {} {} 2>> {}.errlog \n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, modelname,modelname))
        else:
            fid.write('mpiexec -np {} valgrind --leak-check=full  {}/{} {} {}/{} {} 2>> {}.errlog \n'.format(self.nprocs(), self.codepath, executable, solution, self.executionpath, dirname, modelname,modelname))
        if not io_gather: # concatenate the output files
            fid.write('cat {}.outbin.* > {}.outbin'.format(modelname, modelname))
        fid.close()
        fid = open(modelname + '.errlog', 'w') # TODO: Change this to system call (touch <file>)?
        fid.close()
        fid = open(modelname + '.outlog', 'w') # TODO: Change this to system call (touch <file>)?
        fid.close()

    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # since we are running interactively on discover, this doesn't
        # actually upload anything
        subprocess.run(["mkdir","{}".format(dirname)])

        for file in filelist:
            subprocess.run(["cp","./{}".format(file), "{}/{}".format(dirname,file)])
        
        return
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        print('launching solution sequence on discover...')
        if self.runmodel==1:
            # run the model of runmodel=1
            # else, these calls will just write the binary input file
            launchcommand = 'chmod 755 {}.run && ./{}.run 2>> {}.errlog'.format(modelname,modelname,modelname)
            subprocess.call(launchcommand,shell=True)
            subprocess.run(["mv","error.log","./{}/error.log".format(dirname)])
            subprocess.run(["mv","{}.outbin".format(modelname), "./{}/{}.outbin".format(dirname,modelname)])
        # move files to run directory
        subprocess.run(["mv","{}.errlog".format(modelname),"./{}/{}.errlog".format(dirname,modelname)])
        subprocess.run(["mv","{}.outlog".format(modelname), "./{}/{}.outlog".format(dirname,modelname)])
        #subprocess.run(["mv","{}.run".format(modelname), "./{}/{}.run".format(dirname,modelname)])
        return

    def Download(self, dirname, filelist):  # {{{
        # DO NOTHING
        return
    # }}}

def export_discover(md,filename=None,delete_rundir=False):
    # have to change cluster to generic before exporting netCDF because for some reason
    # the .outbin file saves md.cluster as class generic
    # otherwise, loadmodel(filename.nc) throws an error (in loadvars.py module)
    if filename == None:
        filename = './{}/{}.nc'.format(md.private.runtimename,md.miscellaneous.name)
    print('writing model results to netCDF file {}'.format(filename))	
    md.cluster = generic()
    with contextlib.redirect_stderr(io.StringIO()):
    	export_netCDF(md,filename)
    if delete_rundir == True and len(md.private.runtimename)>0:
        subprocess.run(["rm","-rf","./{}".format(md.private.runtimename)])
    return
 
