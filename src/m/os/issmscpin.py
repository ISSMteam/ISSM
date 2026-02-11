import os
import shutil
import subprocess
from MatlabFuncs import *


def issmscpin(host, login, port, path, packages=1):
    """issmscpin get files from host

    Usage:
        issmscpin(host, packages, path)

        bracketstyle:   1 - \\{\\}    (escaped; default)
                        2 - {}      (not escaped)
    """

    # Get hostname
    hostname = oshostname()

    # If hostname and host are the same, do a simple copy
    if strcmpi(hostname, host):
        for package in packages:
            try:
                shutil.copy(os.path.join(path, package), os.getcwd()) # keep going, even if success == 0
            except OSError as e:
                pass
    else:
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
            subproc_cmd = 'scp -P {} {}@localhost:{} {}'.format(port, login, fileliststr, os.getcwd())
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()
            if errs != '':
                # List expansion is a bashism. Try again with '-OT'.
                subproc_cmd = 'scp -OT -P {} {}@localhost:{} {}'.format(port, login, fileliststr, os.getcwd())
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()
        else:
            subproc_cmd = 'scp {}@{}:{} {}'.format(login, host, fileliststr, os.getcwd())
            subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            outs, errs = subproc.communicate()
            if errs != '':
                # List expansion is a bashism. Try again with '-OT'.
                subproc_cmd = 'scp -OT {}@{}:{} {}'.format(login, host, fileliststr, os.getcwd())
                subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                outs, errs = subproc.communicate()

        # Check scp worked
        if errs != '':
            raise OSError('issmscpin error message: {}'.format(errs))
        for package in packages:
            if not os.path.exists(os.path.join('.', package)):
                raise OSError('issmscpin error message: could not scp {}'.format(package))
