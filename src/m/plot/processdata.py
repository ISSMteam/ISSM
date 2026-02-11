import numpy as np
from project2d import project2d
from DepthAverage import DepthAverage

def processdata(md, data, options):
    """PROCESSDATA - process data to be plotted

    datatype = 1 -> elements
    datatype = 2 -> nodes
    datatype = 3 -> node quivers
    datatype = 4 -> P1 patch
    datatype = 5 -> P0 patch
    datatype = 6 -> edges

    Usage:
    data, datatype = processdata(md, data, options)

    See also: PLOTMODEL, PROCESSMESH
    """
    # Initialization and grabbing auxiliaries {{{
    # check format
    if (len(data) == 0 or (len(data) == 1 and not isinstance(data, dict) and np.isnan(data).all())):
        raise ValueError("processdata error message: 'data' provided is empty")
    # get the shape
    if 'numberofvertices2d' in dir(md.mesh):
        numberofvertices2d = md.mesh.numberofvertices2d
    else:
        numberofvertices2d = np.nan

    try:
        numberofedges = md.mesh.numberofedges
    except AttributeError:
        numberofedges = np.nan

    if options.exist('amr'):
        step = options.getfieldvalue('amr', 0)
        numberofvertices = len(md.results.TransientSolution[step].MeshX)
        numberofelements = np.shape(md.results.TransientSolution[step].MeshElements)[0]
    else:
        numberofvertices = md.mesh.numberofvertices
        numberofelements = md.mesh.numberofelements

    procdata = np.copy(data)
    #initialize datatype
    datatype = 0
    # get datasize
    if np.ndim(procdata) == 1:
        datasize = (np.shape(procdata)[0], 1)
    elif np.ndim(procdata) == 2:
        datasize = np.shape(procdata)
    elif np.ndim(procdata) == 3:
        if np.shape(procdata)[0] == 2:
            #treating a dim two list that needs to be stacked
            procdata = np.hstack((procdata[0, :, :], procdata[1, :, :]))
            datasize = np.shape(procdata)
        else:
            raise ValueError('data list contains more than two vectore, we can not cope with that')
    else:
        raise ValueError('data passed to plotmodel has bad dimensions; check that column vectors are rank - 1')
    # }}}

    # log {{{
    if options.exist('log'):
        #cutoff = options.getfieldvalue('log', 1)
        #procdata[np.where(procdata < cutoff)] = cutoff
        #NOTE: Nothing to be done
        pass
    # }}}

    # quiver plot {{{
    if datasize[1] > 1 and datasize[0] != numberofvertices + 1:
        if datasize[0] == numberofvertices and datasize[1] in [2, 3]:
            datatype = 3
            if md.mesh.dimension() == 3:
                if datasize[1] == 2:
                    data = np.hstack(data, np.zeros((datasize[0])))
                elif datasize[1] > 3:
                    raise ValueError('plotmodel error message: data should have two or three columns of length md.mesh.numberofvertices for a quiver plot')
        else:
            #we should have a patch
            print("Assuming that data provided is a patch")
            datatype = 4
            index = md.mesh.elements
            if np.shape(data)[1] < np.shape(index)[1]:
                    raise ValueError('plotmodel error message: data should have more columns than vertices per elements to plot a patch')
            procdata = np.zeros((numberofvertices))
            procdata[md.mesh.elements -1] = data[:, 0:np.shape(index)[1]]
            datasize = [numberofvertices, 1]

    # }}}

    # element data{{{
    if datasize[0] == numberofelements and datasize[1] == 1:
        #initialize datatype if non patch
        if datatype != 4 and datatype != 5:
            datatype = 1
        # AMR {{{
        if options.exist('amr'):
            nonan = np.nonzero(~np.isnan(md.results.TransientSolution[step].MeshElements))
            procdata = procdata[nonan]
        # }}}
        # mask {{{
        if options.exist('mask'):
            flags = options.getfieldvalue('mask')
            hide = np.invert(flags)
            if np.size(flags) == numberofvertices:
                if hasattr(np, 'isin'): #Numpy 2017+
                    tmp = np.isin(index, np.where(hide))
                else: #For backward compatibility
                    tmp = np.in1d(index, np.where(hide))
                EltMask = np.asarray([np.any(tmp) for index in md.mesh.elements - 1])
                procdata = np.ma.array(procdata, mask=EltMask)
                options.addfielddefault('cmap_set_bad', 'w')
            elif np.size(flags) == numberofelements:
                procdata = np.ma.array(procdata, mask=hide)
                options.addfielddefault('cmap_set_bad', 'w')
            else:
                print('Warning: processdata.py: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements')
        # }}}

    # }}}

    # node data {{{
    elif datasize[0] in [numberofvertices, numberofvertices2d] and datasize[1] == 1:
        datatype = 2
        # AMR {{{
        if options.exist('amr'):
            nonan = np.nonzero(~np.isnan(md.results.TransientSolution[step].MeshX))
            procdata = procdata[nonan]
        # }}}
        # Mask {{{
        if options.exist('mask'):
            flags = options.getfieldvalue('mask')
            hide = np.invert(flags)
            if np.size(flags) == numberofvertices:
                procdata = np.ma.array(procdata, mask=hide)
                options.addfielddefault('cmap_set_bad', 'w')
            elif np.size(flags) == numberofelements:
                NodeMask = np.zeros(np.shape(md.mesh.x), dtype=bool)
                HideElt = md.mesh.elements[np.where(hide)[0]] - 1
                NodeMask[HideElt] = True
                procdata = np.ma.array(procdata, mask=NodeMask)
                options.addfielddefault('cmap_set_bad', 'w')
            else:
                print('Warning: processdata.py: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements')
        # }}}
    # }}}

    # edge data {{{
    elif datasize[0] in [numberofedges] and datasize[1] == 1:
        datatype = 6
        procdata = np.zeros((md.mesh.numberofelements))
        repeat = np.zeros((md.mesh.numberofelements))
        for index, edge in enumerate(md.mesh.edges[:, -2:]):
            procdata[edge - 1] += data[index] * np.asarray(edge - 1 > -1, dtype=int)
            repeat[edge - 1] += np.asarray(edge - 1 > -1, dtype=int)
        procdata = procdata / repeat

        # }}}

    # layer procesiong? {{{
    if options.getfieldvalue('layer',0)>=1:
        procdata=project2d(md,data,options.getfieldvalue('layer'))

    if options.getfieldvalue('depthaverage',0):
        raise Exception('Error: "depthaverage" option is not supported in Python.')
        procdata=DepthAverage(md,data) #project onto 2d mesh
    # }}}

    # spc time series {{{
    elif datasize[0] == numberofvertices + 1:
        datatype = 2
        spccol = options.getfieldvalue('spccol', 0)
        print('multiple-column spc field; specify column to plot using option "spccol"')
        print(('column ', spccol, ' plotted for time: ', procdata[-1, spccol]))
        procdata = procdata[0:-1, spccol]

    # }}}

    # convert rank -2 array to rank -1 {{{
    if np.ndim(procdata) == 2 and np.shape(procdata)[1] == 1:
        procdata = procdata.reshape(-1, )
    # }}}

    #  process NaN's if any {{{
    nanfill = options.getfieldvalue('nan', -9999)
    if np.any(np.isnan(procdata)):
        if options.exist('caxis'):
            [lb, ub] = options.getfieldvalue('caxis')
        else:
            lb = np.nanmin(procdata)
            ub = np.nanmax(procdata)
            if lb == ub:
                lb = lb - 0.5
                ub = ub + 0.5
                nanfill = lb - 1
            options.addfielddefault('caxis', [lb, ub])

        procdata[np.isnan(procdata)] = nanfill
        procdata = np.ma.array(procdata, mask=np.isnan(procdata))
        if nanfill < lb:
            options.addfielddefault('cmap_set_under', 'r')
        elif nanfill > ub:
            options.addfielddefault('cmap_set_over', 'k')
        if nanfill < ub and nanfill > lb:
            print(("WARNING: nan's treated as", nanfill, "by default. Which is in your data interval, change it with ['nan', value] in plotmodel options"))
    # }}}

    #  if datatype is still zero, error out {{{
    if datatype == 0:
        raise ValueError("processdata error: data provided not recognized or not supported")
    else:
        return procdata, datatype
    # }}}
