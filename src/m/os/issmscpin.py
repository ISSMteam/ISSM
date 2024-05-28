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
        if ispc():
            #use the putty project pscp.exe: it should be in the path.
            #get ISSM_DIR variable
            if 'ISSM_DIR_WIN' in os.environ:
                ISSM_DIR = os.environ['ISSM_DIR_WIN'][1:-2]
            else:
                raise OSError("issmscpin error message: could not find ISSM_DIR_WIN environment variable.")
            username = eval(input('Username: (quoted string) '))
            key = eval(input('Key: (quoted string) '))
            for package in packages:
                try:
                    subprocess.check_call('%s/externalpackages/ssh/pscp.exe -l "%s" -pw "%s" %s:%s %s' % (ISSM_DIR, username, key, host, os.path.join(path, package), os.getcwd()), shell=True)
                except CalledProcessError as e:
                    raise CalledProcessError("issmscpin error message: could not call putty pscp due to ")

        else:
            #just use standard unix scp
            #string to copy multiple files using scp:
            string = "'{" + ','.join([str(x) for x in packages]) + "}'"
            if port:
                subprocess.call('scp -P {} {}@localhost:{} {}/. '.format(port, login, os.path.join(path, string), os.getcwd()), shell=True)
            else:
                subprocess.call('scp -T {}@{}:{} {}/.'.format(login, host, os.path.join(path, string), os.getcwd()), shell=True)
    #check scp worked
            for package in packages:
                if not os.path.exists(os.path.join('.', package)):
                    raise OSError("issmscpin error message: could not call scp on *nix system for file '{}'".format(package))
