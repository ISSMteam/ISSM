from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import time
import collections
from inspect import isclass
from os import path, remove


class ResTable:
    def __init__(self):
        self.data = []
        self.sizes = []

    def update(self, stepvar):
        #if we have a scalar we just add it to the end
        #we save the size of the current step for further treatment
        if len(np.shape(stepvar)) == 0:
            self.sizes.append(1)
            self.data.append(stepvar)
        # if it is an array we add the values one by one
        #we save the size of the current step for further treatment
        else:
            self.sizes.append([np.shape(stepvar)])
            stackdat = np.squeeze(stepvar.flatten())
            self.data = np.hstack((self.data, stackdat))

    def finalize(self, rows):
        #we have more scalars than steps, so we have an array
        if len(self.data) > rows:
            datasize = np.squeeze(self.sizes)
            maxsize = []
            try:
                dims = np.arange(np.shape(datasize)[1])
                for dim in dims:
                    maxsize.append(np.nanmax(datasize[:, dim]))
            except IndexError:
                if datasize.ndim == 0:
                    maxsize.append(datasize)
                else:
                    maxsize.append(np.nanmax(datasize[:]))
            findim = np.insert(maxsize, 0, rows)
            #first check if all steps are the same size
            if datasize.ndim == 0:
                SameSize = True
            else:
                SameSize = np.sum(np.abs(datasize - datasize[0])) == 0
            if SameSize:
                #same size for all steps, just reshape
                return np.reshape(self.data, newshape=(findim))
            else:
                #different sizes at each steps, first create a table big enough for the biggest step
                startpoint = 0
                outdat = np.nan * np.ones(findim)
                for step in range(rows):
                    #slicer is the data slice in the final array
                    slicer = [slice(0, d) for d in datasize[step, :]]
                    slicer = np.insert(slicer, 0, step)
                    curlen = int(np.prod(datasize[step, :]))
                    outdat[tuple(slicer)] = np.reshape(self.data[startpoint:startpoint + curlen], newshape=(datasize[step, :]))
                    startpoint += curlen
                #outmasked = ma.masked_array(outdat, mask=np.where(np.isnan(outdat), 1, 0))
                return outdat
        #as much scalars as steps (or less) so just one value per step
        else:
            return np.squeeze(np.asarray(self.data))


def export_netCDF(md, filename):  # {{{
    #verbosity of the code, 0 is no messages, 5 is chatty
    verbose = 0
    if path.exists(filename):
        print('File {} allready exist'.format(filename))
        newname = input('Give a new name or "delete" to replace: ')
        if newname == 'delete':
            remove(filename)
        else:
            print(('New file name is {}'.format(newname)))
            filename = newname
    #create file and define it
    NCData = Dataset(filename, 'w', format='NETCDF4')
    NCData.description = 'Results for run' + md.miscellaneous.name
    NCData.history = 'Created ' + time.ctime(time.time())
    # define netCDF dimensions
    #grab time from Transient if it exists
    try:
        StepNum = len(md.results.TransientSolution)
    except AttributeError:  #no transient so just one timestep
        StepNum = 1
    except TypeError:  #this isnot a result so no results in there
        StepNum = 0
    TimeDim = NCData.createDimension('Time', StepNum)  # time is first
    DimDict = {len(TimeDim): 'Time'}
    dimindex = 1
    UnlimDim = NCData.createDimension('Unlim', None)  # unlimited dimension if needed
    DimDict[len(UnlimDim)] = {'Inf'}
    dimindex = 2
    #add mesh related dimension that we know are needed
    dimlist = [2, 40, md.mesh.numberofelements, md.mesh.numberofvertices, np.shape(md.mesh.elements)[1]]
    dimnames = ['DictDummy', 'StringLength', 'EltNum', 'VertNum', 'VertPerElt']
    try:
        dimlist = dimlist + [md.mesh.numberofedges]
        dimnames = dimnames + ['EdgeNum']
    except AttributeError:
        #no edges on this mesh, we fix it at 0
        dimlist += [0]
        dimnames += ['EdgeNum']
    if verbose > 0:
        print('===Creating dimensions ===')
    for i, newdim in enumerate(dimlist):
        if newdim not in list(DimDict.keys()):
            dimindex += 1
            NewDim = NCData.createDimension(dimnames[i], newdim)
            DimDict[len(NewDim)] = dimnames[i]

    typelist = [bool, str, int, float, complex,
                collections.OrderedDict,
                np.int64, np.ndarray, np.float64]

    # get all model classes and create respective groups
    groups = dict.keys(md.__dict__)
    if verbose > 0:
        print('===Creating and populating groups===')
    for group in groups:
        if verbose > 1:
            print('===Now treating {}==='.format(group))
        if group == 'qmu':
            print("qmu is skipped until it is more stable")
            continue
        NCgroup = NCData.createGroup(str(group))
        # In each group gather the fields of the class
        try:
            fields = dict.keys(md.__dict__[group].__dict__)
        except AttributeError:
            print("WARNING: md.{} as no fields, we skip it.".format(group))
            continue
        # looping on fields in each group
        for field in fields:
            Var = md.__dict__[group].__dict__[field]
            # Special treatment for list fields
            if type(Var) == list:
                StdList = False
                if len(Var) == 0:
                    StdList = True  #this is an empty list
                else:
                    #returns False for exotic types (typicaly results)
                    StdList = type(Var[0]) in typelist
                klass = type(md.__dict__[group]).__module__ + '.' + type(md.__dict__[group]).__name__
                NCgroup.__setattr__('classtype', klass)
                if StdList:  # this is a standard or empty list just proceed
                    if verbose > 4:
                        print("=££=creating var for {}.{} with classtype : {}".format(group, field, klass))
                    Var = SqueezeVar(Var)
                    DimDict, ncvar = CreateVar(NCData, Var, field, NCgroup, DimDict)
                    if ncvar is not None:
                        FillVar(ncvar, Var)
                else:  # this is a list of fields, specific treatment needed (usually results or outputdefinitions)
                    if verbose > 4:
                        print("=??=we have a list of fields for {}.{} with classtype : {}".format(group, field, klass))
                    Listsize = len(Var)
                    if group == 'results':  #for results we reshape the datas following time rather than subgrouping
                        Subgroup = NCgroup.createGroup(str(field))
                        try:
                            #take the class of the first element to define nc class and get the list of variables
                            klass = type(md.__dict__[group].__dict__[field][0]).__module__ + '.' + type(md.__dict__[group].__dict__[field][0]).__name__
                            Subgroup.__setattr__('classtype', klass)
                            subfields = dict.keys(md.__dict__[group].__dict__[field][0].__dict__)
                        except (IndexError, AttributeError):
                            klass = type(md.__dict__[group].__dict__[field]).__module__ + '.' + type(md.__dict__[group].__dict__[field]).__name__
                            Subgroup.__setattr__('classtype', klass)
                            subfields = dict.keys(md.__dict__[group].__dict__[field].__getitem__(0))
                        for subfield in subfields:
                            if subfield not in ['errlog', 'outlog']:
                                StackedVar = ResTable()
                                #first loop over the field (result type) to find the index of the last subfield (variable)
                                for listindex in range(0, Listsize):
                                    try:
                                        Var = md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__[subfield]
                                        lastindex = listindex + 1
                                    except AttributeError:
                                        Var = md.__dict__[group].__dict__[field].__getitem__(listindex)[subfield]
                                    except KeyError:
                                        #Some fields only exist for the first step
                                        lastindex = listindex
                                        continue
                                    #Add the  subfield at the current step
                                    Var = SqueezeVar(Var)
                                    StackedVar.update(Var)
                                if verbose > 4:
                                    print("=@@=creating var for {}.{}.{}".format(group, field, subfield))
                                    print("last index of the list is {}".format(lastindex))
                                StackedVar = SqueezeVar(StackedVar.finalize(int(lastindex)))
                                DimDict, ncvar = CreateVar(NCData, StackedVar, subfield, Subgroup, DimDict)
                                #and fill it up
                                if ncvar is not None:
                                    FillVar(ncvar, StackedVar)
                    elif group == 'outputdefinition':  #for outputdefinition we keep a subgroup format
                        for listindex in range(0, Listsize):
                            Subgroupname = str(md.__dict__[group].__dict__[field][listindex].definitionstring)
                            Subgroup = NCgroup.createGroup(Subgroupname)
                            klass = type(md.__dict__[group].__dict__[field][listindex]).__module__ + '.' + type(md.__dict__[group].__dict__[field][listindex]).__name__
                            Subgroup.__setattr__('classtype', klass)
                            subfields = dict.keys(md.__dict__[group].__dict__[field][listindex].__dict__)
                            for subfield in subfields:
                                Var = md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__[subfield]
                                Var = SqueezeVar(Var)
                                if verbose > 4:
                                    print("=--=creating var for {}.{}[{}].{}".format(group, field, listindex, subfield))
                                DimDict, ncvar = CreateVar(NCData, Var, subfield, Subgroup, DimDict)
                                #and fill it up
                                if ncvar is not None:
                                    FillVar(ncvar, Var)
                    else:
                        print("WARNING: unknown treatment for md.{}".format(group))
            # No subgroup, we directly treat the variable
            elif type(md.__dict__[group].__dict__[field]) in typelist or field == 'bamg':
                klass = type(md.__dict__[group]).__module__ + '.' + type(md.__dict__[group]).__name__
                NCgroup.__setattr__('classtype', klass)
                Var = md.__dict__[group].__dict__[field]
                Var = SqueezeVar(Var)
                if verbose > 4:
                    print("====creating var for {}.{}".format(group, field))
                DimDict, ncvar = CreateVar(NCData, Var, field, NCgroup, DimDict)
                if ncvar is not None:
                    FillVar(ncvar, Var)
            # empty field, do nothing
            elif md.__dict__[group].__dict__[field] is None:
                print('field md.{}.{} is None'.format(group, field))
            # if it is a masked array
            elif type(md.__dict__[group].__dict__[field]) is np.ma.core.MaskedArray:
                klass = type(md.__dict__[group]).__module__ + '.' + type(md.__dict__[group]).__name__
                NCgroup.__setattr__('classtype', klass)
                Var = md.__dict__[group].__dict__[field].data
                Var = SqueezeVar(Var)
                if verbose > 4:
                    print("=++=creating var for {}.{}".format(group, field))
                DimDict, ncvar = CreateVar(NCData, Var, field, NCgroup, DimDict)
                if ncvar is not None:
                    FillVar(ncvar, Var)
            # this is an issm class
            elif isclass(type(md.__dict__[group].__dict__[field])):
                if type(md.__dict__[group].__dict__[field]).__name__ == 'solution':
                    #for results we reshape the datas following time rather than subgrouping
                    Listsize = len(md.__dict__[group].__dict__[field])
                    Subgroup = NCgroup.createGroup(str(field))
                    try:
                        #take the class of the first element to define nc class and get the list of variables
                        klass = type(md.__dict__[group].__dict__[field][0]).__module__ + '.' + type(md.__dict__[group].__dict__[field][0]).__name__
                        Subgroup.__setattr__('classtype', klass)
                        subfields = dict.keys(md.__dict__[group].__dict__[field][0].__dict__)
                    except (IndexError, AttributeError):
                        klass = type(md.__dict__[group].__dict__[field]).__module__ + '.' + type(md.__dict__[group].__dict__[field]).__name__
                        Subgroup.__setattr__('classtype', klass)
                        subfields = dict.keys(md.__dict__[group].__dict__[field].__getitem__(0))
                    for subfield in subfields:
                        if subfield not in ['errlog', 'outlog']:
                            StackedVar = ResTable()
                            for listindex in range(0, Listsize):
                                try:
                                    Var = md.__dict__[group].__dict__[field].__getitem__(listindex).__dict__[subfield]
                                    lastindex = listindex + 1
                                except AttributeError:
                                    Var = md.__dict__[group].__dict__[field].__dict__[subfield]
                                    lastindex = listindex
                                except KeyError:
                                    #Some fields only exist for the first step
                                    lastindex = listindex
                                    break
                                Var = SqueezeVar(Var)
                                StackedVar.update(Var)
                            if verbose > 4:
                                print("=$$=creating var for {}.{}.{}".format(group, field, subfield))
                                print("last index of the list is {}".format(lastindex))
                            StackedVar = SqueezeVar(StackedVar.finalize(int(lastindex)))
                            DimDict, ncvar = CreateVar(NCData, StackedVar, subfield, Subgroup, DimDict)
                            #and fill it up
                            if ncvar is not None:
                                FillVar(ncvar, StackedVar)
                elif type(md.__dict__[group].__dict__[field]).__name__ == 'dict':
                    # designed for a dict in dummy but might be used elsewhere
                    # there is no subgroup
                    klass = type(md.__dict__[group]).__module__ + '.' + type(md.__dict__[group]).__name__
                    NCgroup.__setattr__('classtype', klass)
                    Var = md.__dict__[group].__dict__[field]
                    Var = SqueezeVar(Var)
                    if verbose > 4:
                        print("=WW=creating var for {}.{}".format(group, field))
                    DimDict, ncvar = CreateVar(NCData, Var, field, NCgroup, DimDict)
                    if ncvar is not None:
                        FillVar(ncvar, Var)
                else:
                    klass = type(md.__dict__[group]).__module__ + '.' + type(md.__dict__[group]).__name__
                    NCgroup.__setattr__('classtype', klass)
                    Subgroup = NCgroup.createGroup(str(field))
                    klass = type(md.__dict__[group].__dict__[field]).__module__ + '.' + type(md.__dict__[group].__dict__[field]).__name__
                    Subgroup.__setattr__('classtype', klass)
                    subfields = dict.keys(md.__dict__[group].__dict__[field].__dict__)
                    for subfield in subfields:
                        if str(subfield) not in ['errlog', 'outlog']:
                            Var = md.__dict__[group].__dict__[field].__dict__[subfield]
                            Var = SqueezeVar(Var)
                            if verbose > 4:
                                print("+==+creating var for {}.{}.{}".format(group, field, subfield))
                            DimDict, ncvar = CreateVar(NCData, Var, subfield, Subgroup, DimDict)
                            if ncvar is not None:
                                FillVar(ncvar, Var)
            else:
                print("WARNING, md.{}.{} is not treated as it does not fall in one of the existing cases.".format(group, field))

    NCData.close()

# }}}


def CreateVar(NCData, var, field, Group, DimDict, *SupDim):  # {{{
    #=================================================================
    # Define the variables
    #=================================================================
    # grab type
    try:
        val_type = str(var.dtype)
        if val_type.startswith('<U'):
            val_type = 'stringarray'
    except AttributeError:
        val_type = type(var)

    # grab dimension
    if val_type in [collections.OrderedDict, dict]:
        val_shape = len(var.keys())
        val_dim = 2
    else:
        val_shape = np.shape(var)
        val_dim = np.shape(val_shape)[0]
    TypeDict = {float: 'f8',
                'float64': 'f8',
                np.float64: 'f8',
                int: 'i8',
                'int64': 'i8',
                np.int64: 'i8',
                str: str,
                dict: str}

    # Now define and fill up variable
    # treating scalar string or bool as atribute
    if val_type in [str, bool]:
        if field == 'name':  # it looks like netCDF does not like attributes that are called "name"
            field = 'varname'
        Group.__setattr__(str(field), str(var))
        ncvar = None
    # numpy array of strings
    elif val_type == "stringarray":
        #if all strings are the same set it as an attribute
        try:
            samestring = all(var == var[0])
        except IndexError:
            #Only one string
            samestring = True
        if samestring:
            if field == 'name':
                field = 'varname'
            try:
                Group.__setattr__(str(field), str(var[0]))
            except IndexError:
                Group.__setattr__(str(field), str(var))
            ncvar = None
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            ncvar = Group.createVariable(str(field), str, dimensions=dimensions, zlib=True)
    # treating list as string table
    elif val_type == list:
        # try to get the type from the first element
        try:
            nctype = TypeDict[type(var[0])]
        except IndexError:
            nctype = str  # most probably an empty list take str for that

        if val_shape in [(), (0,), 0]:
            ncvar = Group.createVariable(str(field), nctype, zlib=True)
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            ncvar = Group.createVariable(str(field), nctype, dimensions=dimensions)
    # treating  dict as string tables
    elif val_type in [collections.OrderedDict, dict]:
        if val_shape in [(), (0,), 0]:
            ncvar = Group.createVariable(str(field), str, zlib=True)
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            ncvar = Group.createVariable(str(field), str, dimensions=dimensions)
    # treating bool as integers
    elif val_type == 'bool':
        if val_shape in [(), (0,), 0]:
            ncvar = Group.createVariable(str(field), int, zlib=True)
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            ncvar = Group.createVariable(str(field), int, dimensions=dimensions, zlib=True)
    # Now dealing with doubles, we convert them to int if possible
    elif val_type in [float, 'float64', np.float64]:
        try:
            #check if we are integer and under C long overflow also skip empty arrays
            IsInt = np.sum(np.mod(var, 1)) == 0 and np.all(abs(var) < 2147483647) and len(var) > 0
        except TypeError:
            #check if we are integer and under C long overflow
            IsInt = np.mod(var, 1) == 0 and abs(var) < 2147483647
        if IsInt:
            val_type = 'int64'
        if val_shape in [(), (0,), 0] and not SupDim:
            ncvar = Group.createVariable(str(field), TypeDict[val_type], zlib=True)
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            if SupDim:
                dimensions = SupDim + dimensions
            ncvar = Group.createVariable(str(field), TypeDict[val_type], dimensions=dimensions, zlib=True)
    elif val_type in [int, 'int64']:
        if val_shape in [(), (0,), 0] and not SupDim:
            ncvar = Group.createVariable(str(field), TypeDict[val_type], zlib=True)
        else:
            dimensions, DimDict = GetDim(NCData, val_shape, val_type, DimDict, val_dim)
            if SupDim:
                dimensions = SupDim + dimensions
            ncvar = Group.createVariable(str(field), TypeDict[val_type], dimensions=dimensions, zlib=True)
    else:
        print(('WARNING type "{}" is unknown for "{}.{}"'.format(val_type, Group.name, field)))
        ncvar = None
    return DimDict, ncvar
# }}}


def FillVar(ncvar, invar, *UnlimIndex):  # {{{
    #=================================================================
    # Define the variables
    #=================================================================
    # grab type
    try:
        val_type = str(invar.dtype)
        if val_type.startswith('<U'):
            val_type = 'stringarray'
    except AttributeError:
        val_type = type(invar)
    # grab dimension
    if val_type in [collections.OrderedDict, dict]:
        val_shape = len(invar)
    else:
        val_shape = np.shape(invar)

    # Now fill up variable
    # treating list as string table
    if val_type == list:
        if val_shape == 0:
            ncvar = []
        else:
            for elt in range(0, val_shape[0]):
                ncvar[elt] = invar[elt]
    # writing string table
    elif val_type == "stringarray":
        for elt in range(0, val_shape[0]):
            ncvar[elt] = invar[elt]
    # treating bool tables as string tables
    elif val_type == 'bool':
        for elt in range(0, val_shape[0]):
            ncvar[elt] = int(invar[elt])  #str(invar[elt])
    # treating dictionaries as tables of strings
    elif val_type in [collections.OrderedDict, dict]:
        for elt, key in enumerate(dict.keys(invar)):
            ncvar[elt, 0] = key
            ncvar[elt, 1] = str(invar[key])  # converting to str to avoid potential problems
    # Now dealing with numeric variables
    elif val_type in [float, 'float64', np.float64, int, 'int64']:
        try:
            nan_val = np.isnan(invar)
            if nan_val.all():
                naned = 'NaN'
            else:
                naned = invar
        except TypeError:  # type does not accept nan, get vallue of the variable
            naned = invar

        if UnlimIndex:
            if len(val_shape) == 0:
                ncvar[UnlimIndex] = naned
            elif len(val_shape) == 1:
                ncvar[UnlimIndex, :] = naned
            elif len(val_shape) == 2:
                ncvar[UnlimIndex, :, :] = naned
            elif len(val_shape) == 3:
                ncvar[UnlimIndex, :, :, :] = naned
            else:
                print('WARNING: dimension not supported')
        else:
            ncvar[:] = naned
    else:
        print(('WARNING type "{}" is unknown'.format(val_type)))
    return
# }}}


def GetDim(NCData, val_shape, val_type, DimDict, val_dim):  #{{{
    # ============================================================================
    # retriev the dimension tuple from a dictionnary
    # ============================================================================
    output = []
    if val_type in [collections.OrderedDict, dict]:  # dealling with a dictionnary
        try:
            output = [str(DimDict[val_shape])]  # first try to get the coresponding dimension if ti exists
            output = output + [DimDict[2]]  # DictDummy is 2 to treat with dict
        except KeyError:
            index = len(DimDict) + 1  # if the dimension does not exist, increment naming
            NewDim = NCData.createDimension('DimNum' + str(index), val_shape)  # create dimension
            DimDict[len(NewDim)] = 'DimNum' + str(index)  # and update the dimension dictionary
            output = [str(DimDict[val_shape])] + [DimDict[2]]  # now proceed with the shape of the value
    else:
        # loop on dimensions
        for dim in range(0, val_dim):  # loop on the dimensions
            try:
                output = output + [str(DimDict[val_shape[dim]])]  # test if the dimension allready exist
            except KeyError:  # if not create it
                if (val_shape[dim]) > 0:
                    index = len(DimDict) + 1
                    NewDim = NCData.createDimension('DimNum' + str(index), (val_shape[dim]))
                    DimDict[len(NewDim)] = 'DimNum' + str(index)
                    output = output + [str(DimDict[val_shape[dim]])]
    return tuple(output), DimDict
# }}}


def SqueezeVar(Var):  # {{{
    vardim = len(np.shape(Var))
    if vardim > 1:
        Var = np.squeeze(Var)

    return Var
# }}}


def grow(self, row):  # {{{
    np.append(self.data, row)
# }}}
