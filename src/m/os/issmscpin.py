import os
import shutil
import subprocess
from MatlabFuncs import *


def issmscpin(host, login, port, path, packages):
    """ISSMSCPIN get packages from host, using scp on unix, and pscp on windows

    Usage:
        issmscpin(host, packages, path)
    """

    #first get hostname
    hostname = oshostname()
    #first be sure packages are not in the current directory, this could conflict with pscp on windows.
    #remove warnings in case the files do not exist
    for package in packages:
        try:
            os.remove(package)
        except OSError as e:
            pass
    #if hostname and host are the same, do a simple copy
    if strcmpi(hostname, host):  #hostname == host:
        for package in packages:
            try:
                shutil.copy(os.path.join(path, package), os.getcwd())  #keep going, even if success = 0
            except OSError as e:
                pass
    else:
        #just use standard unix scp string to copy multiple files using scp
        filelist = [os.path.join(directory, x) for x in packages]
        fileliststr = ' '.join([str(x) for x in filelist])
        if port:
            subprocess.call('scp -OT -P {} {}@localhost:"{}" {}'.format(port, login, fileliststr, os.getcwd()), shell=True)
        else:
            subprocess.call('scp -OT {}@{}:"{}" {}'.format(login, host, fileliststr, os.getcwd()), shell=True)
        #check scp worked
        for package in packages:
            if not os.path.exists(os.path.join('.', package)):
                raise OSError("issmscpin error message: could not call scp on *nix system for file '{}'".format(package))
