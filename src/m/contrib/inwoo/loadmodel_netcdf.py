#!/usr/bin/env python3
import warnings
from model import *
from pairoptions import pairoptions
import netCDF4
import copy
import pydoc
from results import solution

__all__ = ['loadmodel_netcdf']
def str2bool(v:str): # {{{
    '''
    Okay, convert string type bool, True or False, to bool
    '''
    if not isinstance(v,str):
        raise Exception('ERROR: Given type of input value (={type(v)}) is not string.')

    if v in ['True','1']:
        return True
    else:
        return False
    # }}}

def loadmodel_netcdf(filename, verbose:bool=0):
    '''
    Explain
    -------
    Import model dataset and results from netcdf format. Netcdf format ISSM model result is generated with "export_NetCDF".

    Example
    -------
    >>> import loadmodel_netcdf import loadmodel_netcdf
    >>> md = loadmodel_netcdf('./model.nc') 

    Parameters
    ----------
    verbose:bool (default: 0)
        Show progress with verbosity.

    Returns
    -------
    md: model class in ISSM.
    '''
    # initialize model class..
    md = model();

    # check file name type
    if not filename.endswith('.nc'):
        raise Exception('ERROR: we cannot load %s. This file should be ".nc"(netCDF) file format.')

    # loadmodel
    nc = netCDF4.Dataset(filename, 'r')

    # check groups..
    groupNames = np.array(list(nc.groups.keys()))
    for group, groupName in zip(nc.groups, groupNames):

        # configure variables in class...
        try:
            if (not nc[groupName].variables == {}) | \
                    (groupName in ['private','debug','miscellaneous','verbose','outputdefinition']):
                    #(nc[groupName].classtype == 'results.results') | \
                    #(nc[groupName].classtype == 'private.private') | \
                    #(nc[groupName].classtype == 'debug.debug') | \
                    #(nc[groupName].classtype == 'miscellaneous.miscellaneous') | \
                    #(nc[groupName].classtype == 'verbose.verbose') | \
                    #(nc[groupName].classtype == 'outputdefinition.outputdefinition'):
                mclass = nc[groupName].classtype.split('.')[0]
            elif groupName in ['results']:
                mclass = 'results'
            else: # skip if we cannot configure specific group.
                '''
                Current status
                Some class will be skipeed.
                '''
                warnings.warn(f'WARNING: Skip group = {groupName}')
                continue
        except:
            print(nc.groups)
            raise Exception('ERROR: we cannot find "%s" group.'%(group))

        # get class group.. such as mesh. geometry, inversion, initialization, material etc...
        if hasattr(nc.groups[group],'classgroup'):
            gclass = nc[group].classgroup
        elif hasattr(nc.groups[group],'classtype'):
            gclass = nc.groups[group].classtype
        else:
            gclass = None

        if groupName == 'results': # results special.
            # Loading results with "solution" in ISSM class.
            for solutionname in nc['results'].groups.keys():
                print('loadmodel_netcdf: assign md.results.%s'%(solutionname))
                _solution = solution()
                for fieldname in nc['results'][solutionname].variables.keys():
                    #print(solutionname, fieldname)
                    setattr(_solution, fieldname, np.array(nc['results'][solutionname][fieldname]))

                setattr(md.results, solutionname, _solution)
        elif groupName == 'private': 
            # For private
            md.private.runtimename  = nc[groupName].runtimename
            md.private.isconsistent = str2bool(nc[groupName].isconsistent)
            md.private.solution     = nc[groupName].solution
        elif groupName == 'debug':
            md.debug.valgrind  = str2bool(nc[groupName].valgrind)
            md.debug.gprof     = str2bool(nc[groupName].gprof)
            md.debug.profiling = str2bool(nc[groupName].profiling)
        elif groupName == 'verbose':
            tmp = nc[groupName]
            md.verbose.mprocessor  = str2bool(tmp.mprocessor)
            md.verbose.module      = str2bool(tmp.module)
            try:
                md.verbose.solution    = str2bool(tmp.solution)
            except:
                warnings.warn('Skip loading "md.verbose.solution"')
            md.verbose.solver      = str2bool(tmp.solver)
            md.verbose.convergence = str2bool(tmp.convergence)
            try:
                md.verbose.control     = str2bool(tmp.control)
            except:
                warnings.warn('Skip loading "md.verbose.control"')
            md.verbose.qmu         = str2bool(tmp.qmu)
            md.verbose.autodiff    = str2bool(tmp.autodiff)
            md.verbose.smb         = str2bool(tmp.smb)
        elif groupName == 'miscellaneous':
            # For "miscellaneous"
            md.miscellaneous.notes = nc[groupName].notes
            md.miscellaneous.name  = nc[groupName].varname
        elif groupName == 'outputdefinition':
            definitions = md.outputdefinition.definitions
            if definitions == 'emtpycell':
                md.outputdefinition.definitions = []
            else:
                md.outputdefinition.definitions = definitions
        elif groupName in 'mmemasstransport':
            # NOTE: mmemasstransport is not implemented in python.
            warnings.warn('WARNING: "mmemasstransport" is not implemented in Python.')
            continue
        elif (not mclass in ['verbose',
                                'hpc_simba','fourierlove',
                                ]) & isinstance(gclass,str): # skip specific class..
            # show process..
            if verbose:
                print('    class: {} / classtype: {}'.format(gclass,mclass))

            # get class...
            mod=getattr(__import__(mclass),mclass)()

            # get class
            setattr(md, group, mod)

            # assign loaded calss to md.(field)
            variables = list(nc.groups[group].variables)
            for variable in variables:
                value = nc.groups[group][variable]
                try:
                    if value.shape == ():
                        setattr(getattr(md,group),variable,value[:])
                    else:
                        setattr(getattr(md,group),variable,np.array(value))
                except:
                    value = nc['%s/%s'%(group,variable)]
                    print(value)
                    print(value.shape)
                    raise Exception(f'ERROR: we cannot find {variable}')

    # return outputs.
    return md
