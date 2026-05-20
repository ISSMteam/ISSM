"""Default method implementations shared across HPC cluster classes.

Most remote cluster classes (sherlock, frontera, saga, etc.) share identical
implementations of UploadQueueJob, LaunchQueueJob, and Download.  Rather than
duplicating those blocks in every file, each cluster can delegate to the
corresponding function here.

Usage (inside a cluster class method)::

    import cluster_defaults

    def UploadQueueJob(self, modelname, dirname, filelist):
        cluster_defaults.UploadQueueJob(self, modelname, dirname, filelist)

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):
        cluster_defaults.LaunchQueueJobSbatch(self, modelname, dirname, filelist, restart, batch)

    def Download(self, dirname, filelist):
        cluster_defaults.Download(self, dirname, filelist)

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

    issmscpout(cluster.name, cluster.executionpath, cluster.login, cluster.port, [dirname + '.tar.gz'])
# }}}


def LaunchQueueJobSbatch(cluster, modelname, dirname, filelist, restart, batch):  # {{{
    """Standard SLURM sbatch launch: unpack the tar and submit the queue script via SSH."""
    from helpers import isempty
    if not isempty(restart):
        launchcommand = 'cd {} && cd {} && sbatch {}.queue'.format(
            cluster.executionpath, dirname, modelname)
    else:
        launchcommand = (
            'cd {} && rm -rf ./{} && mkdir {} && cd {} && '
            'mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && sbatch {}.queue'
            .format(cluster.executionpath, dirname, dirname, dirname,
                    dirname, dirname, modelname))
    issmssh(cluster.name, cluster.login, cluster.port, launchcommand)
# }}}


def Download(cluster, dirname, filelist):  # {{{
    """Copy output files from the cluster back to the current directory via scp."""
    directory = '{}/{}/'.format(cluster.executionpath, dirname)
    issmscpin(cluster.name, cluster.login, cluster.port, directory, filelist)
# }}}
