from collections import OrderedDict
# Hack to keep python 2 compatibility
try:
    from dbm import whichdb # Python 3
except ImportError:
    from whichdb import whichdb # Python 2
from re import findall, split
import shelve
from netCDF4 import Dataset, chartostring
import numpy as np
import numpy.ma as ma
from importlib import import_module
from model import *


def loadvars(*args, **kwargs):
    """LOADVARS - function to load variables from a file

    This function loads one or more variables from a file. The names of the
    variables must be supplied. If more than one variable is specified, it may
    be done with a list of names or a dictionary of name as keys. The output
    type will correspond to the input type. All the variables in the file may
    be loaded by specifying only the file name.

    Usage:
        a = loadvars('shelve.dat', 'a')
        [a, b] = loadvars('shelve.dat', ['a', 'b'])
        nvdict = loadvars('shelve.dat', {'a':None, 'b':None})
        nvdict = loadvars('shelve.dat')
    """

    filename = ''
    nvdict = {}
    verbose = 0  # 0 for silent 5 for chatty

    if len(args) >= 1 and isinstance(args[0], str):
        filename = args[0]
        if not filename:
            filename = '/tmp/shelve.dat'
    else:
        raise TypeError("Missing file name.")

    if len(args) >= 2 and isinstance(args[1], str):  # (filename, name)
        for name in args[1:]:
            nvdict[name] = None
    elif len(args) == 2 and isinstance(args[1], list):  # (filename, [names])
        for name in args[1]:
            nvdict[name] = None
    elif len(args) == 2 and isinstance(args[1], dict):  # (filename, {names:values})
        nvdict = args[1]
    elif len(args) == 1:  #  (filename)
        pass
    else:
        raise TypeError("Unrecognized input arguments.")

    timeindex = False
    SteadySols = ['ThermalSolution', 'HydrologySolution', 'StressbalanceSolution']

    for key, value in kwargs.items():
        if key == 'singletime':
            timeindex = value
        if key == 'singleres':
            resname = value

    if whichdb(filename):   #We used python pickle for the save
        print("Loading variables from file {}.".format(filename))
        my_shelf = shelve.open(filename, 'r')  # 'r' for read - only
        if nvdict:
            for name in list(nvdict.keys()):
                try:
                    nvdict[name] = my_shelf[name]
                    print(("Variable '%s' loaded." % name))
                except KeyError:
                    value = None
                    print("Variable '{}' not found.".format(name))

        else:
            for name in list(my_shelf.keys()):
                nvdict[name] = my_shelf[name]
                print(("Variable '%s' loaded." % name))
        my_shelf.close()
    else:  #We used netcdf for the save
        try:
            NCFile = Dataset(filename, mode='r')
            NCFile.close()
        except RuntimeError:
            raise IOError("File '{}' not found.".format(filename))

        classtype, classtree = netCDFread(filename)
        nvdict['md'] = model()
        NCFile = Dataset(filename, mode='r')
        for mod in dict.keys(classtype):
            #==== First we create the model structure  {{{
            if verbose > 0:
                print(' ==== Now treating classtype {}'.format(mod))
            if mod not in classtree.keys():
                print("WARNING: {} classe is not in the model anymore and will be omited.".format(mod))
            elif np.size(classtree[mod]) > 1:
                # this points to a subclass (results.TransientSolution for example)
                curclass = NCFile.groups[classtree[mod][0]].groups[classtree[mod][1]]
                if verbose > 0:
                    print("    ==> {} is of class {}".format(mod, classtype[mod]))
                if classtype[mod][0] == 'results.solutionstep':  #Treating results {{{
                    keylist = [key for key in curclass.groups]
                    #that is the current treatment
                    #here we have a more NC approach with time being a dimension
                    listtype = split(r'\.', classtype[mod][0])[1]
                    try:
                        soltype = str(getattr(curclass, 'SolutionType'))
                    except AttributeError:
                        #might be an older format try that instead :
                        soltype = str(getattr(curclass, 'sOLUTIONtYPE'))
                    if len(NCFile.dimensions['Time']) == 1 or soltype in SteadySols:
                        nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = getattr(classtype[mod][1], listtype)()
                        Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
                    else:
                        #Time dimension is in all the variables so we take that as stepnumber for the results
                        if timeindex:   #we load only the last result to save on time and memory
                            nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = [getattr(classtype[mod][1], listtype)()]
                            Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
                        else:
                            setattr(nvdict['md'].__dict__[classtree[mod][0]], classtree[mod][1], getattr(classtype[mod][1], 'solution')([]))
                            for i in range(max(1, len(NCFile.dimensions['Time']))):
                                nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]].steps.append(getattr(classtype[mod][1], 'solutionstep')())
                            Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]][:]
                # }}}
                elif "results" in mod and classtype[mod][0] == 'list':  #this is the old style of results where every step has a group{{{
                    keylist = [key for key in curclass.groups]
                    #one group per step so use that in place of time
                    stepnum = len(NCFile.groups[classtree[mod][0]].groups[classtree[mod][1]].groups)
                    #we need to redefine classtype from list to result
                    listtype = 'results'
                    classtype[mod].append(__import__(listtype))
                    if stepnum == 1:
                        nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = getattr(classtype[mod][1], listtype)()
                        Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
                    else:
                        if timeindex:   #we load only the last result to save on time and memory
                            nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = [getattr(classtype[mod][1], listtype)()]
                            Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
                        else:
                            setattr(nvdict['md'].__dict__[classtree[mod][0]], classtree[mod][1], getattr(classtype[mod][1], 'solution')([]))
                            for i in range(max(1, stepnum)):
                                nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]].steps.append(getattr(classtype[mod][1], 'solutionstep')())
                            Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]][:]
                    # }}}
                #elif classtype[mod][0] == 'massfluxatgate.massfluxatgate':  #this is for output definitions {{{
                elif mod.startswith('outputdefinition'):  #this is for output definitions {{{
                    defname = split('Output|[0-9]+', classtree[mod][1])[1] + 's'
                    defindex = int(findall('[0-9]+', classtree[mod][1])[0])
                    outdeftype = split(r'\.', classtype[mod][0])[0]
                    nvdict['md'].__dict__[classtree[mod][0]].__dict__[defname].append(getattr(classtype[mod][1], outdeftype)())
                    Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[defname][defindex - 1]
                # }}}
                elif classtype[mod][0] == 'collections.OrderedDict':  #Treating multiple toolkits {{{
                    nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = getattr(classtype[mod][1], 'OrderedDict')
                    Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
                else:
                    if verbose > 0:
                        print("    Using the default for md.{}.{}, is that right??".format(classtree[mod][0], classtree[mod][1]))
                    try:
                        modulename = split(r'\.', classtype[mod][0])[0]
                        if verbose > 0:
                            print("    trying to import {} from {}".format(classtype[mod][0], modulename))
                        nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]] = getattr(classtype[mod][1], modulename)()
                    except AttributeError:
                        print("WARNING: md.{}.{} is not initialized, hopefully that was done in the main group:".format(classtree[mod][0], classtree[mod][1]))
                    Tree = nvdict['md'].__dict__[classtree[mod][0]].__dict__[classtree[mod][1]]
            elif classtype[mod][0] == 'SMBgemb.SMBgemb':
                curclass = NCFile.groups[classtree[mod][0]]
                modulename = split(r'\.', classtype[mod][0])[0]
                nvdict['md'].__dict__[mod] = getattr(classtype[mod][1], modulename)(nvdict['md'].__dict__['mesh'], nvdict['md'].__dict__['geometry'])
                Tree = nvdict['md'].__dict__[classtree[mod][0]]
            else:
                curclass = NCFile.groups[classtree[mod][0]]
                modulename = split(r'\.', classtype[mod][0])[0]
                nvdict['md'].__dict__[mod] = getattr(classtype[mod][1], modulename)()
                Tree = nvdict['md'].__dict__[classtree[mod][0]]
            if verbose > 0:
                print("    for {} Tree is a {} with len {}".format(mod, Tree.__class__.__name__, len(curclass.groups)))
            # }}}
            #==== Then we populate it {{{
            #for i in range(0, max(1, len(curclass.groups))):
            if len(curclass.groups) > 0:  #that is presumably only for old style NC where each result step had its own group
                if timeindex:
                    if timeindex < 0:
                        groupclass = [curclass.groups[keylist[len(curclass.groups) - timeindex]]]
                    else:
                        groupclass = [curclass.groups[keylist[timeindex]]]
                else:
                    groupclass = [curclass.groups[key] for key in keylist]
            else:
                groupclass = [curclass]
            for groupindex, listclass in enumerate(groupclass):
                try:
                    soltype = str(getattr(listclass, 'SolutionType'))
                except AttributeError:
                    soltype = 'NoSol'
                #==== We deal with Variables {{{
                for var in listclass.variables:
                    if not resname or var == resname:
                        if var not in ['errlog', 'outlog']:
                            varval = listclass.variables[str(var)]
                            vardim = varval.ndim
                            if verbose > 0:
                                print("    ==> treating var {} of dimension {}".format(var, vardim))
                            #There is a special treatment for results to account for its specific structure
                            #that is the new export version where time is a named dimension
                            NewFormat = 'Time' in NCFile.dimensions
                            if type(Tree) == list:  # and NewFormat:
                                if timeindex:
                                    if NewFormat:
                                        if vardim == 0:
                                            try:
                                                Tree[0].__dict__[str(var)] = varval[timeindex].data
                                            except IndexError:
                                                print('WARNING: No data on index {} for {} reverting to last time.'.format(timeindex, str(var)))
                                                Tree[0].__dict__[str(var)] = varval[-1].data
                                        elif vardim == 1:
                                            try:
                                                Tree[0].__dict__[str(var)] = varval[timeindex].data
                                            except IndexError:
                                                print('WARNING: No data on index {} for {} reverting to last time.'.format(timeindex, str(var)))
                                                Tree[0].__dict__[str(var)] = varval[-1].data
                                        elif vardim == 2:
                                            Tree[0].__dict__[str(var)] = varval[timeindex, :].data
                                        elif vardim == 3:
                                            Tree[0].__dict__[str(var)] = varval[timeindex, :, :].data
                                        else:
                                            print('table dimension greater than 3 not implemented yet')
                                    elif soltype in SteadySols:
                                        Tree.__dict__[str(var)] = varval[:].data
                                    else:  #old format had step sorted in difeerent group so last group is last time
                                        Tree[0].__dict__[str(var)] = varval[:].data
                                else:
                                    if NewFormat:
                                        incomplete = 'Time' not in varval.dimensions and soltype not in SteadySols
                                        if incomplete:
                                            try:
                                                chosendim = varval.dimensions[0]
                                                timelist = np.arange(0, len(NCFile.dimensions[chosendim]))
                                                print('WARNING, {} is not present on every times, we chose {}({}) as the dimension to write it with'.format(var, chosendim, len(NCFile.dimensions[chosendim])))
                                            except IndexError:
                                                #just one step, so no dimension, we just put it on the first solutionstep
                                                timelist = [0]
                                        elif soltype in SteadySols:
                                            timelist = [0]
                                        else:
                                            timelist = np.arange(0, len(NCFile.dimensions['Time']))
                                        if soltype in SteadySols:
                                            Tree.__dict__[str(var)] = varval[:].data
                                        else:
                                            for t in timelist:
                                                if verbose > 5:
                                                    print("filing step {} for {}".format(t, var))
                                                if vardim == 0:
                                                    Tree[t].__dict__[str(var)] = varval[:].data
                                                elif vardim == 1:
                                                    stepval = ma.masked_array(varval[t].data, mask=np.where(np.isnan(varval[t]), 1, 0))
                                                    Tree[t].__dict__[str(var)] = ma.compressed(stepval)
                                                elif vardim == 2:
                                                    stepval = ma.masked_array(varval[t, :].data, mask=np.where(np.isnan(varval[t, :]), 1, 0))
                                                    Tree[t].__dict__[str(var)] = ma.compressed(stepval)
                                                elif vardim == 3:
                                                    stepval = ma.masked_array(varval[t, :, :].data, mask=np.where(np.isnan(varval[t, :, :]), 1, 0))
                                                    Tree[t].__dict__[str(var)] = ma.compressed(stepval).reshape((stepval.count(0)[0], stepval.count(1)[0]))
                                                else:
                                                    print('table dimension greater than 3 not implemented yet')
                                    else:
                                        if verbose > 0:
                                            print("filing step {} for {}".format(groupindex, var))
                                        Tree[groupindex].__dict__[str(var)] = varval[:].data
                            else:
                                if vardim == 0:  #that is a scalar
                                    if str(varval[0]) in ['', '--', 'emptycell']:  #no value
                                        Tree.__dict__[str(var)] = []
                                    elif varval[0] == 'True':  #treatin bool
                                        Tree.__dict__[str(var)] = True
                                    elif varval[0] == 'False':  #treatin bool
                                        Tree.__dict__[str(var)] = False
                                    else:
                                        Tree.__dict__[str(var)] = varval[0].item()

                                elif vardim == 1:  #that is a vector
                                    if verbose > 0:
                                        print("   for variable {} type is {}".format(str(var), varval.dtype))
                                    if varval.dtype == str:
                                        if varval.shape[0] == 1:
                                            Tree.__dict__[str(var)] = [str(varval[0]), ]
                                        elif 'True' in varval[:] or 'False' in varval[:]:
                                            Tree.__dict__[str(var)] = np.asarray([V == 'True' for V in varval[:]], dtype=bool)
                                        else:
                                            Tree.__dict__[str(var)] = [str(vallue) for vallue in varval[:]]
                                    elif varval.dtype == "|S1":  #that is for matlab chararcter arrays
                                        stringlist = chartostring(varval[:])
                                        Tree.__dict__[str(var)] = [stringlist.tolist(), ]
                                    else:
                                        try:
                                            #some thing specifically require a list
                                            mdtype = type(Tree.__dict__[str(var)])
                                        except KeyError:
                                            mdtype = float
                                        if mdtype == list:
                                            Tree.__dict__[str(var)] = [mdval for mdval in varval[:]]
                                        else:
                                            Tree.__dict__[str(var)] = varval[:].data

                                elif vardim == 2:
                                    #dealling with dict
                                    if verbose > 0:
                                        print("   for variable {} type is {}".format(str(var), varval.dtype))
                                    if varval.dtype == str:  #that is for dictionaries
                                        if any(varval[:, 0] == 'toolkit'):  #toolkit definition have to be first
                                            Tree.__dict__[str(var)] = OrderedDict([('toolkit', str(varval[np.where(varval[:, 0] == 'toolkit')[0][0], 1]))])
                                            strings1 = [str(arg[0]) for arg in varval if arg[0] != 'toolkits']
                                            strings2 = [str(arg[1]) for arg in varval if arg[0] != 'toolkits']
                                            Tree.__dict__[str(var)].update(list(zip(strings1, strings2)))
                                        else:
                                            strings1 = [str(arg[0]) for arg in varval]
                                            strings2 = [str(arg[1]) for arg in varval]
                                            Tree.__dict__[str(var)] = OrderedDict(list(zip(strings1, strings2)))
                                    elif varval.dtype == "|S1":  #that is for matlab chararcter arrays
                                        stringlist = chartostring(varval[:, :])
                                        stringlist = [string.strip() for string in stringlist]
                                        Tree.__dict__[str(var)] = stringlist
                                    else:
                                        if type(Tree) == list:
                                            t = indexlist[i]
                                            if listtype == 'dict':
                                                Tree[t][str(var)] = varval[:, :].data
                                            else:
                                                Tree[t].__dict__[str(var)] = varval[:, :].data
                                        else:
                                            Tree.__dict__[str(var)] = varval[:, :].data
                                elif vardim == 3:
                                    if varval.dtype == "|S1":  #that is for matlab chararcter arrays
                                        #most likely that is a toolkit dictionar so should be treated as such
                                        #first we convert the character table to strings
                                        stringtable = []
                                        for i in range(np.shape(varval)[0]):
                                            stringtable.append([chartostring(varval[i, 0, :]), chartostring(varval[i, 1, :])])
                                        stringtable = np.asarray(stringtable, dtype=str)
                                        Tree.__dict__[str(var)] = OrderedDict([('toolkit', str(varval[np.where(stringtable[:, 0] == 'toolkit')[0][0], 1]))])
                                        strings1 = [str(arg[0]) for arg in stringtable if arg[0] != 'toolkits']
                                        strings2 = [str(arg[1]) for arg in stringtable if arg[0] != 'toolkits']
                                        Tree.__dict__[str(var)].update(list(zip(strings1, strings2)))
                                    else:
                                        Tree.__dict__[str(var)] = varval[:, :, :].data
                                else:
                                    print('table dimension greater than 3 not implemented yet')
                # }}}
                #==== And with atribute {{{
                for attr in listclass.ncattrs():
                    if verbose > 0:
                        print("      ==> treating attribute {}".format(attr))
                    if attr != 'classtype':  #classtype is for treatment, don't get it back
                        if attr == 'varname':
                            attribute = 'name'
                        else:
                            attribute = attr
                        if type(Tree) == list:
                            if verbose > 0:
                                print("        printing with index 0")
                            if listtype == 'dict':
                                Tree[0][attribute] = str(listclass.getncattr(attr))
                            else:
                                Tree[0].__dict__[attribute] = str(listclass.getncattr(attr))
                        else:
                            if listclass.getncattr(attr) == 'True':
                                Tree.__dict__[attribute] = True
                            elif listclass.getncattr(attr) == 'False':
                                Tree.__dict__[attribute] = False
                            elif listclass.getncattr(attr) == 'emptycell':
                                Tree.__dict__[attribute] = []
                            else:
                                Tree.__dict__[attribute] = str(listclass.getncattr(attr))
                # }}}
            # }}}
        NCFile.close()
    if len(args) >= 2 and isinstance(args[1], str):  # (value)
        value = [nvdict[name] for name in args[1:]]
        return value

    elif len(args) == 2 and isinstance(args[1], list):  # ([values])
        value = [nvdict[name] for name in args[1]]
        return value

    elif (len(args) == 2 and isinstance(args[1], dict)) or (len(args) == 1):  # ({names:values})
        return nvdict


def netCDFread(filename):
    print(('Opening {} for reading '.format(filename)))
    NCData = Dataset(filename, 'r')
    class_dict = {}
    class_tree = {}

    for group in NCData.groups:
        if len(NCData.groups[group].groups) > 0:
            for subgroup in NCData.groups[group].groups:
                classe = str(group) + '.' + str(subgroup)
                grpclass = str(getattr(NCData.groups[group].groups[subgroup], 'classtype'))
                class_dict[classe] = [grpclass, ]
                if class_dict[classe][0] not in ['dict', 'list', 'cell']:
                    try:
                        modulename = split(r'\.', class_dict[classe][0])[0]
                        class_dict[classe].append(import_module(modulename))
                    except ModuleNotFoundError:
                        #submodule probably has a different name
                        modulename = str(getattr(NCData.groups[group].groups[subgroup], 'classtype'))
                        print("WARNING importing {} rather than {}".format(modulename, class_dict[classe][0]))
                        class_dict[classe].append(import_module(modulename))
                class_tree[classe] = [group, subgroup]
        else:
            classe = str(group)
            try:
                class_dict[classe] = [str(getattr(NCData.groups[group], 'classtype')), ]
                if class_dict[classe][0] not in ['dict', 'list', 'cell']:
                    modulename = split(r'\.', class_dict[classe][0])[0]
                    try:
                        class_dict[classe].append(import_module(modulename))
                        class_tree[classe] = [group, ]
                    except ModuleNotFoundError:
                        print("WARNING: module {} does not exist anymore and is skipped".format(modulename))

            except AttributeError:
                print(('group {} is empty'.format(group)))
    NCData.close()
    return class_dict, class_tree
