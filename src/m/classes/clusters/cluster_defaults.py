"""Default method implementations shared across HPC cluster classes.

Most remote cluster classes (sherlock, frontera, saga, etc.) share identical
implementations of UploadQueueJob, LaunchQueueJob, and Download.  Rather than
duplicating those blocks in every file, each cluster can delegate to the
corresponding function here.

Clusters whose methods differ from these defaults (e.g. bbftp transfers,
interactive path variants, qsub/PBS schedulers, or local execution) should
keep their own implementations unchanged.
"""

import os
import subprocess

from issmdir import issmdir
from issmscpin import issmscpin
from issmscpout import issmscpout
from issmssh import issmssh

def UploadQueueJob(cluster, modelname, dirname, filelist):  # {{{
    """Compress the execution directory into a tar.gz and scp it to the cluster.

    filelist contains full paths; ``-C`` is used so only basenames are stored
    in the archive.
    """
    root = os.path.join(issmdir(), 'execution', dirname)
    compressstring = 'tar -C {} -zcf {}.tar.gz'.format(root, dirname)
    for filepath in filelist:
        if not os.path.isfile(filepath):
            raise Exception('File {} not found'.format(filepath))
        compressstring += ' {}'.format(os.path.basename(filepath))
    subprocess.call(compressstring, shell=True)

    port = getattr(cluster, 'port', 0)
    issmscpout(cluster.name, cluster.executionpath, cluster.login, port, [dirname + '.tar.gz'])
# }}}

def LaunchQueueJobSbatch(cluster, modelname, dirname, filelist, restart, batch, fmt):  # {{{
    """Launch a queued job on the cluster via SSH.

    fmt values:
        1 - no scheduler, source .queue directly
        2 - SLURM, use sbatch
        3 - PBS/Torque, use qsub
    """
    from helpers import isempty
    if fmt == 1:
        # No scheduler: source .queue directly
        if not isempty(restart):
            launchcommand = 'source {}/environment.sh && cd {} && cd {} && source {}.queue'.format(
                cluster.etcpath, cluster.executionpath, dirname, modelname)
        else:
            if batch:
                launchcommand = ('source {}/environment.sh && cd {} && rm -rf ./{} && mkdir {} && '
                                 'cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz'
                                 .format(cluster.etcpath, cluster.executionpath, dirname, dirname,
                                         dirname, dirname, dirname))
            else:
                launchcommand = 'source {}/environment.sh && source {}/{}/{}.queue'.format(
                    cluster.etcpath, cluster.executionpath, dirname, modelname)
    elif fmt == 2:
        # SLURM sbatch
        if not isempty(restart):
            launchcommand = 'cd {} && cd {} && sbatch {}.queue'.format(
                cluster.executionpath, dirname, modelname)
        else:
            launchcommand = (
                'cd {} && rm -rf ./{} && mkdir {} && cd {} && '
                'mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && sbatch {}.queue'
                .format(cluster.executionpath, dirname, dirname, dirname,
                        dirname, dirname, modelname))
    elif fmt == 3:
        # PBS/Torque qsub
        if not isempty(restart):
            launchcommand = 'cd {} && cd {} && /PBS/bin/qsub {}.queue'.format(
                cluster.executionpath, dirname, modelname)
        else:
            launchcommand = (
                'cd {} && rm -rf ./{} && mkdir {} && cd {} && '
                'mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && /PBS/bin/qsub {}.queue'
                .format(cluster.executionpath, dirname, dirname, dirname,
                        dirname, dirname, modelname))
    else:
        raise ValueError('fmt={} not supported'.format(fmt))
    port = getattr(cluster, 'port', 0)
    issmssh(cluster.name, cluster.login, port, launchcommand)
# }}}

def Download(cluster, dirname, filelist):  # {{{
    """Copy output files from the cluster back to the current directory via scp."""
    directory = '{}/{}/'.format(cluster.executionpath, dirname)
    port = getattr(cluster, 'port', 0)
    issmscpin(cluster.name, cluster.login, port, directory, filelist)
# }}}
