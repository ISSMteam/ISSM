from datetime import datetime
import os

import numpy as np

from loadresultsfromcluster import loadresultsfromcluster
from pairoptions import pairoptions
from waitonlock import waitonlock


def solveslm(slm, solutionstringi, *args):
    """solveslm - apply solution sequence for this sealevel model

    Usage:
        slm = solveslm(slm,solutionstring,varargin)
        where varargin is a lit of paired arguments of string OR enums

    solution types available comprise:
        - 'Transient'

    extra options:

    Examples:
        slm=solveslm(slm,'Transient');
    """

    # Recover and process solve options
    if solutionstringi.lower() == 'tr' or solutionstringi.lower() == 'transient':
        solutionstring = 'TransientSolution'
    else:
        raise RuntimeError('solutionstring {} not supported!'.format(solutionstringi))

    # Default settings for debugging
    valgrind = 0
    #slm.cluster.interactive = 0
    #valgrind = 1

    # Check consistency
    slm.checkconsistency(solutionstring)

    # Process options
    options = pairoptions('solutionstring', solutionstring, *args)

    # Make sure we request sum of cluster processors
    totalnp = 0
    for i in range(len(slm.icecaps)):
        totalnp = totalnp + slm.icecaps[i].cluster.np
    totalnp = totalnp + slm.earth.cluster.np
    if totalnp != slm.cluster.np:
        raise RuntimeError('sum of all icecaps and earth cluster processors requests should be equal to slm.cluster.np')

    # Recover some fields
    slm.private.solution = solutionstring
    cluster = slm.cluster
    batch = 0
    # Now, go through icecaps, glaciers and earth, and upload all the data independently
    print('solving ice caps first')
    for i in range(len(slm.icecaps)):
        slm.icecaps[i] = solve(slm.icecaps[i], solutionastringi,'batch','yes')
    print('solving earth now')
    slm.earth = solve(slm.earth, solutionstringi, 'batch', 'yes')

    # First, build a runtime name that is unique
    c = datetime.now()
    md.private.runtimename = "%s-%02i-%02i-%04i-%02i-%02i-%02i-%i" % (md.miscellaneous.name, c.month, c.day, c.year, c.hour, c.minute, c.second, os.getpid())

    # Write all input files
    privateruntimenames = []
    miscellaneousnames = []
    nps = []
    for i in range(len(slm.icecaps)):
        privateruntimenames.append(slm.icecaps[i],private.runtimename)
        miscellaneousnames.append(slm.earth.miscellaneous.name)
        nps.append(slm.earth.cluster.np)

    BuildQueueScriptMultipleModels(cluster, slm.private.runtimename, slm.miscellaneous.name, slm.private.solution, privateruntimenames, miscellaneousnames, nps)

    # Upload all required files, given that each individual solution for icecaps and earth model already did
    filelist = [slm.miscellaneous.name + '.queue']
    UploadQueueJob(cluster, slm.miscellaneous.name, slm.private.runtimename, filelist)

    # Launch queue job
    print('launching solution sequence')
    LaunchQueueJob(cluster, slm.miscellaneous.name, slm.private.runtimename, filelist, '', batch)

    # Wait on lock
    if slm.settings.waitonlock > 0:
        islock = waitonlock(slm)
        if islock == 0:  # no results to be loaded
            print('The results must be loaded manually with md = loadresultsfromcluster(md).')
        else: # load results
            if slm.verbose.solution:
                print('loading results from cluster')
            slm = loadresultsfromcluster(slm)

    return slm
