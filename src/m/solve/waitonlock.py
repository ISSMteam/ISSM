import subprocess
import sys
import time
from generic import *
from generic_static import *
# from localpfe import *
from MatlabFuncs import *
from pfe import *

def waitonlock(md):
    """WAITONLOCK - wait for a file

    This routine will return when a file named 'lockfilename' is written to
    disk. Also check for outlog file because it might be written several
    seconds after the lock file.

    If the time limit given in input is exceeded, return 0.

    Usage:
        flag = waitonlock(md)

    TODO:
    - Uncomment import of localpfe and check on cluster type once localpfe.py
    has been translated from localpfe.m.
    """

    # Get lockfilename (lock file) and options
    executionpath = md.cluster.executionpath
    timelimit = md.settings.waitonlock
    cluster = md.cluster

    if isa(cluster, pfe) and cluster.interactive > 1:
        lockfilename = '{}/Interactive{}/{}.lock'.format(executionpath, cluster.interactive, md.miscellaneous.name)
        logfilename = '{}/Interactive{}/{}.outlog'.format(executionpath, cluster.interactive, md.miscellaneous.name)
    # elif isa(cluster, localpfe):
    #     lockfilename = '{}/{}.lock'.format(executionpath, md.miscellaneous.name)
    #     logfilename = '{}/{}.outlog'.format(executionpath, md.miscellaneous.name)
    else:
        lockfilename = '{}/{}/{}.lock'.format(executionpath, md.private.runtimename, md.miscellaneous.name)
        logfilename = '{}/{}/{}.outlog'.format(executionpath, md.private.runtimename, md.miscellaneous.name)

    # If we are using the generic cluster in interactive mode, job is already complete
    if (isa(cluster, generic) and cluster.interactive) or (isa(cluster, generic_static)):
        # We are in interactive mode, no need to check for job completion
        return 1

    # Initialize time and file presence test flag
    elapsedtime = 0
    ispresent = 0
    starttime = time.time()
    print('waiting for {} hold on... (Ctrl+C to exit)'.format(lockfilename))

    # Prepare command if the job is not running on the local machine
    if not strcmpi(oshostname(), cluster.name):
        if cluster.name == 'cloud':
            command = '[ -f {} ] && [ -f {} ] 2>/dev/null'.format(lockfilename, logfilename)
            command = '{} sshmaster {} --user {} \'{}\''.format(starcluster(), cluster.name, cluster.login, command)
        else:
            command = 'ssh -l {}'.format(cluster.login)
            if isprop(cluster, 'idfile') and cluster.idfile != '':
                command += ' -i {}'.format(cluster.idfile)
            if isprop(cluster, 'port') and cluster.port:
                command += ' -p {} localhost'.format(cluster.port);
            else:
                command += ' {}'.format(cluster.name)
            command += ' "[ -f {} ] && [ -f {} ]" 2>/dev/null'.format(lockfilename, logfilename)

    while not ispresent and elapsedtime < timelimit:
        if strcmpi(oshostname(), cluster.name):
            pause(1)
            ispresent = (isfile(lockfilename) and isfile(logfilename))
            elapsedtime = etime(time.time(), starttime) / 60
        else:
            pause(5)
            elapsedtime = etime(time.time(), starttime)
            sys.stdout.write('\rchecking for job completion (time: {} min {} sec)      '.format(floor(elapsedtime / 60), floor(rem(elapsedtime, 60)))) # TODO: After Python 2 is deprecated, we can change this call to print([...], end='')
            elapsedtime = elapsedtime / 60 # Converts time from sec to min
            subproc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()  # NOTE: Need to consume output before checking return code

            if errs != '':
               raise Exception('waitonlock: check for existence of files failed: {}'.format(errs))
            ispresent = not subproc.returncode
            if ispresent:
                print('')

    # Build output
    if elapsedtime > timelimit:
        print('Time limit exceeded. Increase md.settings.waitonlock')
        print('The results must be loaded manually with md = loadresultsfromcluster(md).')
        raise RuntimeError('waitonlock error message: time limit exceeded.')

    return ispresent
