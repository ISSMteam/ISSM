from datetime import datetime
import os

from ismodelselfconsistent import ismodelselfconsistent
from loadresultsfromcluster import loadresultsfromcluster
from marshall import marshall
from pairoptions import pairoptions
from preqmu import *
from waitonlock import waitonlock


def solve(md, solutionstring, *args):
    """solve - apply solution sequence for this model

    Usage:
        md = solve(md, solutionstring, varargin)

    where varargin is a list of paired arguments of string OR enums

    Solution types available comprise:
    - 'Stressbalance'        or 'sb'
    - 'Masstransport'        or 'mt'
    - 'Oceantransport'       or 'oceant'
    - 'Thermal'              or 'th'
    - 'Steadystate'          or 'ss'
    - 'Transient'            or 'tr'
    - 'Balancethickness'     or 'mc'
    - 'BalancethicknessSoft' or 'mcsoft'
    - 'Balancevelocity'      or 'bv'
    - 'BedSlope'             or 'bsl'
    - 'SurfaceSlope'         or 'ssl'
    - 'Hydrology'            or 'hy'
    - 'DamageEvolution'      or 'da'
    - 'Gia'                  or 'gia'
    - 'Love'                 or 'lv'
    - 'Esa'                  or 'esa'
    - 'Sampling'             or 'smp'

    Extra options:
    - loadonly         : do not solve, only load results
    - runtimename      : true or false (default is true); makes name unique
    - checkconsistency : 'yes' or 'no' (default is 'yes'); checks consistency 
                         of model
    - restart          : directory name (relative to the execution directory) 
                         where the restart file is located

    Examples:
        md = solve(md, 'Stressbalance')
        md = solve(md, 'sb')
    """

    # Recover and process solve options
    if solutionstring.lower() == 'sb' or solutionstring.lower() == 'stressbalance':
        solutionstring = 'StressbalanceSolution'
    elif solutionstring.lower() == 'mt' or solutionstring.lower() == 'masstransport':
        solutionstring = 'MasstransportSolution'
    elif solutionstring.lower() == 'oceant' or solutionstring.lower() == 'oceantransport':
        solutionstring = 'OceantransportSolution'
    elif solutionstring.lower() == 'th' or solutionstring.lower() == 'thermal':
        solutionstring = 'ThermalSolution'
    elif solutionstring.lower() == 'st' or solutionstring.lower() == 'steadystate':
        solutionstring = 'SteadystateSolution'
    elif solutionstring.lower() == 'tr' or solutionstring.lower() == 'transient':
        solutionstring = 'TransientSolution'
    elif solutionstring.lower() == 'mc' or solutionstring.lower() == 'balancethickness':
        solutionstring = 'BalancethicknessSolution'
    elif solutionstring.lower() == 'bv' or solutionstring.lower() == 'balancevelocity':
        solutionstring = 'BalancevelocitySolution'
    elif solutionstring.lower() == 'bsl' or solutionstring.lower() == 'bedslope':
        solutionstring = 'BedSlopeSolution'
    elif solutionstring.lower() == 'ssl' or solutionstring.lower() == 'surfaceslope':
        solutionstring = 'SurfaceSlopeSolution'
    elif solutionstring.lower() == 'hy' or solutionstring.lower() == 'hydrology':
        solutionstring = 'HydrologySolution'
    elif solutionstring.lower() == 'da' or solutionstring.lower() == 'damageevolution':
        solutionstring = 'DamageEvolutionSolution'
    elif solutionstring.lower() == 'gia' or solutionstring.lower() == 'gia':
        solutionstring = 'GiaSolution'
    elif solutionstring.lower() == 'lv' or solutionstring.lower() == 'love':
        solutionstring = 'LoveSolution'
    elif solutionstring.lower() == 'esa':
        solutionstring = 'EsaSolution'
    elif solutionstring.lower() == 'smp' or solutionstring.lower() == 'sampling':
        solutionstring = 'SamplingSolution'
    else:
        raise ValueError('solutionstring {} not supported!'.format(solutionstring))
    options = pairoptions('solutionstring', solutionstring, *args)

    # Recover some fields
    md.private.solution = solutionstring
    cluster = md.cluster
    if options.getfieldvalue('batch', 'no') == 'yes':
        batch = 1
    else:
        batch = 0

    # Check model consistency
    if options.getfieldvalue('checkconsistency', 'yes') == 'yes':
        if md.verbose.solution:
            print('checking model consistency')
        ismodelselfconsistent(md)

    # If we are restarting, actually use the provided runtime name
    restart = options.getfieldvalue('restart', '')
    # First, build a runtime name that is unique
    if restart == 1:
        pass # Leave the runtimename as is
    else:
        if not isempty(restart):
            md.private.runtimename = restart
        else:
            if options.getfieldvalue('runtimename', True):
                c = datetime.now()
                md.private.runtimename = '%s-%02i-%02i-%04i-%02i-%02i-%02i-%i' % (md.miscellaneous.name, c.month, c.day, c.year, c.hour, c.minute, c.second, os.getpid())
            else:
                md.private.runtimename = md.miscellaneous.name

    # If running QMU analysis, some preprocessing of Dakota files using model 
    # fields needs to be carried out
    if md.qmu.isdakota:
        md = preqmu(md, options)

    # Do we load results only?
    if options.getfieldvalue('loadonly', False):
        md = loadresultsfromcluster(md)
        return md

    # Write all input files
    marshall(md) # bin file
    md.toolkits.ToolkitsFile(md.miscellaneous.name + '.toolkits') # toolkits file
    cluster.BuildQueueScript(md.private.runtimename, md.miscellaneous.name, md.private.solution, md.settings.io_gather, md.debug.valgrind, md.debug.gprof, md.qmu.isdakota, md.transient.isoceancoupling) # queue file

    # Upload all required files
    modelname = md.miscellaneous.name
    filelist = [modelname + '.bin', modelname + '.toolkits']

    if ispc():
        filelist.append(modelname + '.bat')
    else:
        filelist.append(modelname + '.queue')

    if md.qmu.isdakota:
        filelist.append(modelname + '.qmu.in')

    if isempty(restart):
        print('uploading input files')
        cluster.UploadQueueJob(md.miscellaneous.name, md.private.runtimename, filelist)

    # Launch job
    print('launching solution sequence')
    cluster.LaunchQueueJob(md.miscellaneous.name, md.private.runtimename, filelist, restart, batch)

    # Return if batch
    if batch:
        if md.verbose.solution:
            print('batch mode requested: not launching job interactively')
            print('launch solution sequence on remote cluster by hand')
        return md

    # Wait on lock
    if md.settings.waitonlock > 0:
        # Wait for done file
        done = waitonlock(md)
        if md.verbose.solution:
            print('loading results from cluster')
        md = loadresultsfromcluster(md)
    elif md.settings.waitonlock == 0:
        print('Model results must be loaded manually with md = loadresultsfromcluster(md).')

    return md
