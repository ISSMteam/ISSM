import os
import shutil
import subprocess
from MatlabFuncs import *


def issmscpin(host, login, port, path, packages):
    """issmscpin get packages from host using scp

    Usage:
        issmscpin(host, packages, path)

    NOTE: If users again have issues with file list (i.e.

        {<FILE1>,<FILE2>,...<FILEN>}

    ), note that this a bash'ism and default shell should be checked. View file 
    history for potential fix (i.e. some combination of -O and -T options).
    """

    # First get hostname
    hostname = oshostname()

    # If hostname and host are the same, do a simple copy
    if strcmpi(hostname, host):
        for package in packages:
            try:
                shutil.copy(os.path.join(path, package), os.getcwd()) # keep going, even if success == 0
            except OSError as e:
                pass
    else:
        filelist = [os.path.join(directory, x) for x in packages]
        fileliststr = ' '.join([str(x) for x in filelist])
        if port:
            subprocess.call('scp -P {} {}@localhost:"{}" {}'.format(port, login, fileliststr, os.getcwd()), shell=True)
        else:
            subprocess.call('scp {}@{}:"{}" {}'.format(login, host, fileliststr, os.getcwd()), shell=True)
        # Check scp worked
        for package in packages:
            if not os.path.exists(os.path.join('.', package)):
                raise OSError('issmscpin error message: could not scp {}'.format(package))
