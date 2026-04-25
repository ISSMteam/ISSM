import numpy as np
import os
from re import findall, match, split
from pairoptions import pairoptions
from operator import attrgetter
import MatlabFuncs as m


def checkfield(md, *args):
    """CHECKFIELD - check field consistency

    Used to check model consistency
    Requires:
        'field' or 'fieldname' option. If 'fieldname' is provided, it will
        retrieve it from the model md. (md.(fieldname))
        If 'field' is provided, it will assume the argument following 'field'
        is a numeric array.

    Available options:
        - NaN: 1 if check that there is no NaN
        - size: [lines cols], NaN for non checked dimensions, or 'universal'
        for any input type (nodal, element, time series, etc)
        - > :  greater than provided value
        - >= : greater or equal to provided value
        - < :  smallerthan provided value
        - <=: smaller or equal to provided value
        - < vec:  smallerthan provided values on each vertex
        - timeseries: 1 if check time series consistency (size and time)
        - values: list of acceptable values
        - numel: list of acceptable number of elements
        - cell: 1 if check that is cell
        - empty: 1 if check that non empty
        - message: overloaded error message

    Usage:
        md = checkfield(md, fieldname, options)
    """

    #get options
    options = pairoptions(*args)

    #get field from model
    field = None
    if options.exist('field'):
        field = options.getfieldvalue('field')
        fieldname = options.getfieldvalue('fieldname', 'no fieldname')
    else:
        fieldname = options.getfieldvalue('fieldname')
        fieldnametokens = split(r'\.', fieldname)
        for i in range(len(fieldnametokens)):
            if match(r'(.*?)\[(.*?)\]\[(.*?)\]', fieldnametokens[i]):
                # Handle 2D list and list of dict addressing
                listname = split(r'\[(.*?)\]', fieldnametokens[i])[0]
                listindices = findall(r'\[(.*?)\]', fieldnametokens[i])
                if field == None:
                    field = attrgetter(listname)(md)
                else:
                    field = attrgetter(listname)(field)
                try:
                    listindex = int(listindices[0])
                    field = field[listindex]
                    try:
                        listindex = int(listindices[1])
                    except:
                        listindex = listindices[1][1:-1] # Strip enclosing quote
                    field = field[listindex]
                except:
                    listindex = listindices[0][1:-1] # Strip enclosing quote
                    field = field[listindex]
                    try:
                        listindex = int(listindices[1])
                    except:
                        listindex = listindices[1][1:-1] # Strip enclosing quote
                    field = field[listindex]
            elif match(r'(.*?)\[(.*?)\]', fieldnametokens[i]):
                # Handle list and dict addressing
                listname = split(r'\[(.*?)\]', fieldnametokens[i])[0]
                listindex = findall(r'\[(.*?)\]', fieldnametokens[i])[0]
                if field == None:
                    field = attrgetter(listname)(md)
                else:
                    field = attrgetter(listname)(field)
                try:
                    listindex = int(listindex)
                    field = field[listindex]
                except:
                    field = attrgetter(listindex)(field)
            else:
                if field == None:
                    field = attrgetter(fieldnametokens[i])(md)
                else:
                    field = attrgetter(fieldnametokens[i])(field)

    if isinstance(field, (bool, int, float)):
        field = np.array([field])

    #check empty
    if options.exist('empty'):
        if not field:
            md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is empty'.format(fieldname)))

    #Check size
    if options.exist('size'):
        fieldsize = options.getfieldvalue('size')
        if type(fieldsize) == str:
            if fieldsize == 'universal':
                #Check that vector size will not be confusing for ModelProcessorx
                if (md.mesh.numberofvertices == md.mesh.numberofelements):
                    raise Exception('number of vertices is the same as number of elements')
                elif (md.mesh.numberofvertices + 1 == md.mesh.numberofelements):
                    raise Exception('number of vertices + 1 is the same as number of elements')
                elif (md.mesh.numberofvertices == md.mesh.numberofelements + 1):
                    raise Exception('number of vertices is the same as number of elements + 1')

                #Uniform field
                if (field.shape[0] == 1):
                    if (np.ndim(field) > 1 and np.shape(field)[1] != 1):
                        md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is not supported'.format(fieldname)))

                #vertex oriented input, only one column allowed
                elif (np.shape(field)[0] == md.mesh.numberofvertices):
                    if (np.ndim(field) > 1 and np.shape(field)[1] != 1):
                        md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is not supported'.format(fieldname)))

                #element oriented input, one or more column (patch) is ok
                elif (np.shape(field)[0] == md.mesh.numberofelements):
                    #nothing to do here (either constant per element, or defined on nodes)
                    pass

                #vertex time series
                elif (np.shape(field)[0] == md.mesh.numberofvertices + 1):
                    if (np.ndim(field) > 1 and np.shape(field)[1] <= 1):
                        md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is not supported'.format(fieldname)))

                #element time series
                elif (np.shape(field)[0] == md.mesh.numberofelements + 1):
                    if (np.ndim(field) > 1 and np.shape(field)[1] <= 1):
                        md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is not supported'.format(fieldname)))

                #else not supported
                else:
                    md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' is not supported'.format(fieldname)))

            else:
                raise Exception('fieldsize \'{}\' not supported yet'.format(fieldsize))

        else:
            if len(np.shape(field)) < len(fieldsize):
                if fieldsize[-1] > 1:
                    md = md.checkmessage(options.getfieldvalue('message', 'field {} has size {} but should be size {}'.format(fieldname, np.shape(field), fieldsize)))
                else:
                    #The last dimension is one that follows MATLAB 2D array regulation but usually not what we do in Python, we allow the difference in shape only if the number of element is equal
                    if np.prod(np.shape(field)) != np.prod(fieldsize):
                        md = md.checkmessage(options.getfieldvalue('message', 'field {} has size {} but should be size {}'.format(fieldname, np.shape(field), fieldsize)))
            else:
                for i in range(np.size(fieldsize)):
                    if (not np.isnan(fieldsize[i])) and (np.shape(field)[i] != fieldsize[i]):
                        md = md.checkmessage(options.getfieldvalue('message', 'field {} dimension # {} should be of size {}'.format(fieldname, i, fieldsize[i])))

    #Check numel
    if options.exist('numel'):
        fieldnumel = options.getfieldvalue('numel')
        if (type(fieldnumel) == int and np.size(field) != fieldnumel) or (type(fieldnumel) == list and np.size(field) not in fieldnumel):
            if len(fieldnumel) == 1:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' size should be {}'.format(fieldname, fieldnumel)))
            elif len(fieldnumel) == 2:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' size should be {} or {}'.format(fieldname, fieldnumel[0], fieldnumel[1])))
            else:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' size should be {}'.format(fieldname, fieldnumel)))

    #check NaN
    if options.getfieldvalue('NaN', 0):
        if np.any(np.isnan(field)):
            md = md.checkmessage(options.getfieldvalue('message', 'NaN values found in field \'{}\''.format(fieldname)))

    #check Inf
    if options.getfieldvalue('Inf', 0):
        if np.any(np.isinf(field)):
            md = md.checkmessage(options.getfieldvalue('message', 'Inf values found in field \'{}\''.format(fieldname)))

    #check cell
    if options.getfieldvalue('cell', 0):
        if not isinstance(field, (tuple, list, dict)):
            md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' should be a tuple, list, or dict'.format(fieldname)))

    #check values
    if options.exist('values'):
        fieldvalues = options.getfieldvalue('values')
        if False in m.ismember(field, fieldvalues):
            if len(fieldvalues) == 1:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' value should be \'{}\''.format(fieldname, fieldvalues[0])))
            elif len(fieldvalues) == 2:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' values should be \'{}\' or \'{}\''.format(fieldname, fieldvalues[0], fieldvalues[1])))
            else:
                md = md.checkmessage(options.getfieldvalue('message', 'field \'{}\' should have values in {}'.format(fieldname, fieldvalues)))

    #check greater
    if options.exist('>='):
        lowerbound = options.getfieldvalue('>=')
        if type(lowerbound) is str:
            lowerbound = attrgetter(lowerbound)(md)
        if np.size(lowerbound) > 1:  #checking elementwise
            if any(field < lowerbound):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values above {}'.format(fieldname, lowerbound)))
        else:
            minval = np.nanmin(field)
            if options.getfieldvalue('timeseries', 0):
                minval = np.nanmin(field[:-1])
            elif options.getfieldvalue('singletimeseries', 0):
                if np.size(field) == 1:  #some singletimeseries are just one value
                    minval = field
                else:
                    minval = np.nanmin(field[0])
            elif options.getfieldvalue('mappedtimeseries', 0):
                minval = np.nanmin(field[:-1])

            if minval < lowerbound:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values above {}'.format(fieldname, lowerbound)))

    if options.exist('>'):
        lowerbound = options.getfieldvalue('>')
        if type(lowerbound) is str:
            lowerbound = attrgetter(lowerbound)(md)
        if np.size(lowerbound) > 1:  #checking elementwise
            if any(field <= lowerbound):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values above {}'.format(fieldname, lowerbound)))
        else:
            minval = np.nanmin(field)
            if options.getfieldvalue('timeseries', 0):
                minval = np.nanmin(field[:-1])
            elif options.getfieldvalue('singletimeseries', 0):
                if np.size(field) == 1:  #some singletimeseries are just one value
                    minval = field
                else:
                    minval = np.nanmin(field[0])
            elif options.getfieldvalue('mappedtimeseries', 0):
                minval = np.nanmin(field[:-1])

            if minval <= lowerbound:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values above {}'.format(fieldname, lowerbound)))

    #check smaller
    if options.exist('<='):
        upperbound = options.getfieldvalue('<=')
        if type(upperbound) is str:
            upperbound = attrgetter(upperbound)(md)
        if np.size(upperbound) > 1:  #checking elementwise
            if any(field > upperbound):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values below {}'.format(fieldname, upperbound)))
        else:
            maxval = np.nanmax(field)
            if options.getfieldvalue('timeseries', 0):
                maxval = np.nanmax(field[:-1])
            elif options.getfieldvalue('singletimeseries', 0):
                if np.size(field) == 1:  #some singletimeseries are just one value
                    maxval = field
                else:
                    maxval = np.nanmax(field[0])
            elif options.getfieldvalue('mappedtimeseries', 0):
                maxval = np.nanmax(field[:-1])
            elif hasattr(field, 'fov_forward_indices'):
                maxval = field.fov_forward_indices[0]
            if maxval > upperbound:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values below {}'.format(fieldname, upperbound)))

    if options.exist('<'):
        upperbound = options.getfieldvalue('<')
        if type(upperbound) is str:
            upperbound = attrgetter(upperbound)(md)
        if np.size(upperbound) > 1:  #checking elementwise
            if any(field >= upperbound):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values below {}'.format(fieldname, upperbound)))

        else:
            maxval = np.nanmax(field)
            if options.getfieldvalue('timeseries', 0):
                maxval = np.nanmax(field[:-1])
            elif options.getfieldvalue('singletimeseries', 0):
                if np.size(field) == 1:  #some singletimeseries are just one value
                    maxval = field.copy()
                else:
                    maxval = np.nanmax(field[0])
            elif options.getfieldvalue('mappedtimeseries', 0):
                maxval = np.nanmax(field[:-1])

            if maxval >= upperbound:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have values below {}'.format(fieldname, upperbound)))

    #check file
    if options.getfieldvalue('file', 0):
        if not os.path.exists(field):
            md = md.checkmessage('file provided in {}: {} does not exist'.format(fieldname, field))

    #Check row of strings
    if options.exist('stringrow'):
        if not isinstance(field, list):
            md = md.checkmessage(options.getfieldvalue('message', 'field {} should be a list'.format(fieldname)))

    #Check forcings (size and times)
    if options.getfieldvalue('timeseries', 0):
        if field.shape[0] == md.mesh.numberofvertices or field.shape[0] == md.mesh.numberofelements:
            if np.ndim(field) > 1 and not np.size(field, 1) == 1:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should have only one column as there are md.mesh.numberofvertices lines'.format(fieldname)))
        elif field.shape[0] == md.mesh.numberofvertices + 1 or field.shape[0] == md.mesh.numberofelements + 1:
            if np.ndim(field) > 1 and not all(field[-1, :] == np.sort(field[-1, :])):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} columns should be sorted chronologically'.format(fieldname)))
            if np.ndim(field) > 1 and any(field[-1, 0:-1] == field[-1, 1:]):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} columns must not contain duplicate timesteps'.format(fieldname)))
        else:
            md = md.checkmessage(options.getfieldvalue('message', 'field {} should have md.mesh.numberofvertices or md.mesh.numberofvertices + 1 lines'.format(fieldname)))

    #Check single value forcings (size and times)
    if options.getfieldvalue('singletimeseries', 0):
        if field.shape[0] == 2:
            if not all(field[-1, :] == np.sort(field[-1, :])):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} columns should be sorted chronologically'.format(fieldname)))
            if any(field[-1, 0:-1] == field[-1, 1:]):
                md = md.checkmessage(options.getfieldvalue('message', 'field {} columns must not contain duplicate timesteps'.format(fieldname)))
        elif field.shape[0] == 1:
            if np.ndim(field) > 1 and not np.size(field, 1) == 1:
                md = md.checkmessage(options.getfieldvalue('message', 'field {} should be either a scalar or have 2 lines'.format(fieldname)))
        else:
            md = md.checkmessage(options.getfieldvalue('message', 'field {} should have 2 lines or be a scalar'.format(fieldname)))

    #Check forcings (size and times)
    if options.getfieldvalue('mappedtimeseries', 0):
        if np.ndim(field) > 1 and not all(field[-1, :] == np.sort(field[-1, :])):
            md = md.checkmessage(options.getfieldvalue('message', 'field {} columns should be sorted chronologically'.format(fieldname)))
        if np.ndim(field) > 1 and any(field[-1, 0:-1] == field[-1, 1:]):
            md = md.checkmessage(options.getfieldvalue('message', 'field {} columns must not contain duplicate timesteps'.format(fieldname)))

    return md
