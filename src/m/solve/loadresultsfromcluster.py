import os
import subprocess

from helpers import *
from loadresultsfromdisk import loadresultsfromdisk
from MatlabFuncs import *
from pairoptions import pairoptions


def remove(filename):  #{{{
    try:
        os.remove(filename)
    except OSError:
        print(('WARNING: ' + filename + ' does not exist'))
# }}}


def loadresultsfromcluster(md, *args):  #{{{
    """loadresultsfromcluster - load results of solution sequence from cluster

    Usage:
        md = loadresultsfromcluster(md)
        md = loadresultsfromcluster(md, 'runtimename', runtimename)

        Options include: 'runtimename', 'nolog'

    Example:
        md = loadresultsfromcluster(md, 'runtimename', 'test101-06-15-2021-13-24-18-4883')
    """

    # Process options
    options = pairoptions(*args)
    nolog = options.getfieldvalue('nolog', 0)
    md.private.runtimename = options.getfieldvalue('runtimename', md.private.runtimename)

    # Retrieve cluster, to be able to call its methods
    cluster = md.cluster

    # Download outputs from the cluster
    if not nolog:
        filelist = [md.miscellaneous.name + '.outlog', md.miscellaneous.name + '.errlog']
    else:
        filelist = []
    if md.qmu.isdakota:
        filelist.append(md.miscellaneous.name + '.qmu.err')
        filelist.append(md.miscellaneous.name + '.qmu.out')
        if 'tabular_graphics_data' in fieldnames(md.qmu.params):
            if md.qmu.params.tabular_graphics_data:
                filelist.append('dakota_tabular.dat')
        if md.qmu.output and md.qmu.statistics.method[0]['name'] == 'None':
            if md.qmu.method.method == 'nond_sampling':
                for i in range(md.qmu.method.params.samples):
                    filelist.append(md.miscellaneous.name + '.outbin.' + str(i + 1))
        if md.qmu.statistics.method[0]['name'] != 'None':
            filelist.append(md.miscellaneous.name + '.stats')
    else:
        filelist.append(md.miscellaneous.name + '.outbin')
    cluster.Download(md.private.runtimename, filelist)

    # If we are here, no errors in the solution sequence, call loadresultsfromdisk
    md = loadresultsfromdisk(md, md.miscellaneous.name + '.outbin')

    # Erase the log and output files
    for i in range(len(filelist)):
        filename = filelist[i]
        remove(filename)
    if not ispc():
        remove(md.private.runtimename + '.tar.gz')

    # Erase input file if run was carried out on same platform
    hostname = oshostname()
    if hostname == cluster.name:
        remove(md.miscellaneous.name + '.bin')
        remove(md.miscellaneous.name + '.toolkits')
        if md.qmu.isdakota:
            remove(md.miscellaneous.name + '.qmu.in')
        if not ispc():
            remove(md.miscellaneous.name + '.queue')
        else:
            remove(md.miscellaneous.name + '.bat')

    return md
# }}}
