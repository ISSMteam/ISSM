import os
import subprocess
from MatlabFuncs import *


def issmscpout(host, path, login, port, packages):
    """ISSMSCPOUT send packages to a host, using scp on unix, and pscp on windows

    Usage:
        issmscpout(host, path, packages)
    """

    #get hostname
    hostname = oshostname()

    #if hostname and host are the same, do a simple copy

    if strcmpi(host, hostname):  #host == hostname:
        for package in packages:
            here = os.getcwd()
            os.chdir(path)
            try:
                os.remove(package)
            except OSError:
                pass
            subprocess.call('ln -s %s %s' % (os.path.join(here, package), path), shell=True)
            os.chdir(here)
    else:
        if ispc():
            #use the putty project pscp.exe: it should be in the path.
            #get ISSM_DIR variable
            if 'ISSM_DIR_WIN' in os.environ:
                ISSM_DIR = os.environ['ISSM_DIR_WIN'][1:-2]
            else:
                raise OSError("issmscpout error message: could not find ISSM_DIR_WIN environment variable.")

            username = eval(input('Username: (quoted string) '))
            key = eval(input('Key: (quoted string) '))

            for package in packages:
                try:
                    subprocess.check_call('%s/externalpackages/ssh/pscp.exe -l "%s" -pw "%s" %s %s:%s' % (ISSM_DIR, username, key, package, host, path), shell=True)
                except CalledProcessError as e:
                    raise CalledProcessError("issmscpout error message: could not call putty pscp.")

        else:
            #just use standard unix scp
            #create string of packages being sent
            string = ''
            for package in packages:
                string += ' ' + package
            string += ' '

            if port:
                subprocess.call('scp -P %d %s %s@localhost:%s' % (port, string, login, path), shell=True)
            else:
                subprocess.call('scp %s %s@%s:%s' % (string, login, host, path), shell=True)
