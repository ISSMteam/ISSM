import os
import subprocess
from MatlabFuncs import *


def issmscpout(host, path, login, port, packages):
    """issmscpout send files to host

    Usage:
        issmscpout(host, path, packages)
    """

    # Get hostname
    hostname = oshostname()

    # If hostname and host are the same, do a simple copy
    if strcmpi(host, hostname):
        for package in packages:
            here = os.getcwd()
            os.chdir(path)
            try:
                os.remove(package)
            except OSError:
                pass
            subprocess.call('ln -s %s %s' % (os.path.join(here, package), path), shell=True)
            os.chdir(here)

    #General case, this is not a local machine
    else:
        filelist = [os.path.join(directory, x) for x in packages]
        fileliststr = ' '.join([str(x) for x in filelist])
        if port:
            subproc_cmd = 'scp -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path)
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()
            if errs != '':
                subproc_cmd = 'scp -OT -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path)
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()
        else:
            subproc_cmd = 'scp {} {}@{}:{}'.format(fileliststr, login, host, path)
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()
            if errs != '':
                subproc_cmd = 'scp -OT {} {}@{}:{}'.format(fileliststr, login, host, path)
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()

        # Check scp worked
        if errs != '':
            raise OSError('issmscpin error message: {}'.format(errs))
