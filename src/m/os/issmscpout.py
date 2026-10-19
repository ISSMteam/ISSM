import os
import shlex
import subprocess
import warnings
from MatlabFuncs import *
from issmscpin import haversync, runcmd


def rsyncout(host, login, port, path, packages):  # {{{
    """Push packages with rsync, passing every package as a source.
    """

    if port:
        remote = '{}@localhost:'.format(login)
        sshopt = ' -e {}'.format(shlex.quote('ssh -p {}'.format(port)))
    else:
        remote = '{}@{}:'.format(login, host)
        sshopt = ''

    sources = ' '.join([shlex.quote(package) for package in packages])
    return runcmd('rsync -a{} {} {}'.format(sshopt, sources, remote + shlex.quote(path)))
# }}}


def scpout(host, login, port, path, packages, noescape):  # {{{
    """Push packages with scp, relying on the local shell to expand \\{\\}
    """

    if len(packages) == 1:
        fileliststr = packages[0]
    else:
        fileliststr = r'\{'
        fileliststr += ','.join([package for package in packages])
        fileliststr += r'\}'

        # Remove backslashes if noescape is set
        if noescape:
            fileliststr = fileliststr[1:-2] + fileliststr[-1]

    if port:
        status, errs = runcmd('scp -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path))
        if status != 0:
            # List expansion is a bashism. Try again with '-OT'.
            status, errs = runcmd('scp -OT -P {} {} {}@localhost:{}'.format(port, fileliststr, login, path))
    else:
        status, errs = runcmd('scp {} {}@{}:{}'.format(fileliststr, login, host, path))
        if status != 0:
            # List expansion is a bashism. Try again with '-OT'.
            status, errs = runcmd('scp -OT {} {}@{}:{}'.format(fileliststr, login, host, path))

    return status, errs
# }}}


def issmscpout(host, path, login, port, packages, transfermethod='auto'):
    """issmscpout send files to host

    Usage:
        issmscpout(host, path, login, port, packages, transfermethod)

        transfermethod: 'auto'         - use rsync, fall back to scp if the remote host
                                         does not have rsync installed (default)
                        'rsync'        - use rsync only
                        'scp'          - use scp only, with escaped brackets \\{\\}
                        'scp-noescape' - use scp only, with unescaped brackets {}
                                         (for hosts whose shell must not see backslashes)
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
    usescp = transfermethod.lower() in ('scp', 'scp-noescape')
    if not usescp and not strcmpi(transfermethod, 'rsync') and not strcmpi(transfermethod, 'auto'):
        raise ValueError('issmscpout error message: unknown transfermethod \'{}\''.format(transfermethod))

    # No point in trying rsync if this machine does not have it
    if not usescp and not haversync():
        if strcmpi(transfermethod, 'rsync'):
            raise OSError('issmscpout error message: rsync requested but not found on this machine')
        usescp = True

    if not usescp:
        status, errs = rsyncout(host, login, port, path, packages)

        # rsync returns 127 when the remote host does not have rsync installed.
        # Any other non-zero status is a genuine transfer error and must not
        # trigger a second connection (it would prompt for credentials twice).
        if status == 127 and strcmpi(transfermethod, 'auto'):
            warnings.warn('issmscpout: rsync not found on {}, falling back to scp'.format(host))
            usescp = True

    if usescp:
        status, errs = scpout(host, login, port, path, packages, strcmpi(transfermethod, 'scp-noescape'))

    # Check that the transfer worked
    if status != 0:
        raise OSError('issmscpout error message: {}'.format(errs))
