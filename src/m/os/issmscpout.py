import os
import subprocess
from MatlabFuncs import *


def issmscpout(host, path, login, port, packages, bracketstyle=1):
    """issmscpout send files to host

    Usage:
        issmscpout(host, path, login, port, packages, bracketstyle)

        bracketstyle:   1 - \\{\\}    (escaped; default)
                        2 - {}      (not escaped)
    """

    # Get hostname
    hostname = oshostname()

    # If hostname and host are the same, do a simple copy or symlinks
    if strcmpi(host, hostname):
        here = os.getcwd()
        for package in packages:
            try:
                os.remove(os.path.join(path, package))
            except OSError:
                pass
            subprocess.call('ln -s {} {}'.format(os.path.join(here, package), path), shell=True)
        return 0

    # General case: this is not a local machine
    if len(packages) == 1:
        fileliststr = packages[0]
    else:
        fileliststr = r'\{'
        fileliststr += ','.join([package for package in packages])
        fileliststr += r'\}'

        # Remove backslashes if bracketstyle is 2
        if bracketstyle == 2:
            fileliststr = fileliststr[1:-2] + fileliststr[-1]
    if port:
        subproc_cmd = 'scp -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path)
        subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        outs, errs = subproc.communicate()
        if errs != '':
            # List expansion is a bashism. Try again with '-OT'.
            subproc_cmd = 'scp -OT -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path)
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()
    else:
        subproc_cmd = 'scp {} {}@{}:{}'.format(fileliststr, login, host, path)
        subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        outs, errs = subproc.communicate()
        if errs != '':
            # List expansion is a bashism. Try again with '-OT'.
            subproc_cmd = 'scp -OT {} {}@{}:{}'.format(fileliststr, login, host, path)
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()

    # Check scp worked
    if errs != '':
        raise OSError('issmscpin error message: {}'.format(errs))
