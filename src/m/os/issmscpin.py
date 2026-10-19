import os
import shlex
import shutil
import subprocess
import warnings
from MatlabFuncs import *


_have_rsync = None


def haversync():  # {{{
    """Is rsync available on this machine? (probed once per session)
    """

    global _have_rsync

    if _have_rsync is None:
        _have_rsync = shutil.which('rsync') is not None
    return _have_rsync
# }}}


def runcmd(subproc_cmd):  # {{{
    """Run a transfer command, returning its exit status and stderr.

    The exit status is what decides success: scp and rsync both write progress
    and ssh banners to stderr even when the transfer worked.
    """

    subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = subproc.communicate()
    return subproc.returncode, errs
# }}}


def rsyncin(host, login, port, path, filelist):  # {{{
    """Pull filelist with rsync, one remote source per file.

    All sources share the same host, so rsync uses a single connection. Each
    source is quoted individually, so no list expansion is needed.
    """

    if port:
        remote = '{}@localhost:'.format(login)
        sshopt = ' -e {}'.format(shlex.quote('ssh -p {}'.format(port)))
    else:
        remote = '{}@{}:'.format(login, host)
        sshopt = ''

    sources = ' '.join([remote + shlex.quote(path + '/' + file) for file in filelist])
    return runcmd('rsync -a{} {} {}'.format(sshopt, sources, shlex.quote(os.getcwd())))
# }}}


def scpin(host, login, port, path, filelist, noescape):  # {{{
    """Pull filelist with scp, relying on the remote shell to expand \\{\\}
    """

    if len(filelist) == 1:
        fileliststr = filelist[0]
    else:
        fileliststr = r'\{'
        fileliststr += ','.join([file for file in filelist])
        fileliststr += r'\}'

        # Remove backslashes if noescape is set
        if noescape:
            fileliststr = fileliststr[1:-2] + fileliststr[-1]

    if port:
        status, errs = runcmd('scp -P {} {}@localhost:{}/{} {}'.format(port, login, path, fileliststr, os.getcwd()))
        if status != 0:
            # List expansion is a bashism. Try again with '-OT'.
            status, errs = runcmd('scp -OT -P {} {}@localhost:{}/{} {}'.format(port, login, path, fileliststr, os.getcwd()))
    else:
        status, errs = runcmd('scp {}@{}:{}/{} {}'.format(login, host, path, fileliststr, os.getcwd()))
        if status != 0:
            # List expansion is a bashism. Try again with '-OT'.
            status, errs = runcmd('scp -OT {}@{}:{}/{} {}'.format(login, host, path, fileliststr, os.getcwd()))

    return status, errs
# }}}


def issmscpin(host, login, port, path, filelist=1, transfermethod='auto'):
    """issmscpin get files from host

    Usage:
        issmscpin(host, login, port, path, filelist, transfermethod)

        transfermethod: 'auto'         - use rsync, fall back to scp if the remote host
                                         does not have rsync installed (default)
                        'rsync'        - use rsync only
                        'scp'          - use scp only, with escaped brackets \\{\\}
                        'scp-noescape' - use scp only, with unescaped brackets {}
                                         (for hosts whose shell must not see backslashes)
    """

    # Get hostname
    hostname = oshostname()

    # If hostname and host are the same, do a simple copy
    if strcmpi(hostname, host):
        for file in filelist:
            try:
                shutil.copy(os.path.join(path, file), os.getcwd()) # keep going, even if success == 0
            except OSError as e:
                warnings.warn('Could not copy ' + os.path.join(path, file))
                pass
        return

    usescp = transfermethod.lower() in ('scp', 'scp-noescape')
    if not usescp and not strcmpi(transfermethod, 'rsync') and not strcmpi(transfermethod, 'auto'):
        raise ValueError('unknown transfermethod \'{}\''.format(transfermethod))

    # No point in trying rsync if this machine does not have it
    if not usescp and not haversync():
        if strcmpi(transfermethod, 'rsync'):
            raise OSError('rsync requested but not found on this machine')
        usescp = True

    if not usescp:
        status, errs = rsyncin(host, login, port, path, filelist)

        # rsync returns 127 when the remote host does not have rsync installed.
        # Any other non-zero status is a genuine transfer error and must not
        # trigger a second connection (it would prompt for credentials twice).
        if status == 127 and strcmpi(transfermethod, 'auto'):
            warnings.warn('issmscpin: rsync not found on {}, falling back to scp'.format(host))
            usescp = True

    if usescp:
        status, errs = scpin(host, login, port, path, filelist, strcmpi(transfermethod, 'scp-noescape'))

    # Check that the transfer worked
    if status != 0:
        raise OSError('{}'.format(errs))
    for file in filelist:
        if not os.path.exists(os.path.join('.', file)):
            raise OSError('could not transfer {}'.format(file))
