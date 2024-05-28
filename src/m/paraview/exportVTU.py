import numpy as np
from base64 import b64encode
from os import path, remove, mkdir
from glob import glob


def exportVTU(filename, md, *args, enveloppe=False, fmtout="binary", **kwargs):
    '''
    vtu export
    function exportVTU(filename, md)
    Exports resluts in XML based vtu format for visualisation in Paraview.
    Hopefully it is based on the treatment for export VTK and only the output part is modified.
    (only work for triangle and wedges based on their number of nodes)

    Usage:
    exportVTU('FileName', md)
    exportVTU('FileName', md, 'geometry', 'mesh')
    exportVTU('FileName', md, 'geometry', 'mesh', enveloppe = True)

    DirName is the name of the output directory, each timestep then has it
    own file ('Timestep.vtkX.vtk') with X the number of the output step
    enveloppe is an option keeping only the enveloppe of the md (it is False by default)

    Options:
        - clipping : allows to reduce your domain (cliping=[Xmin, Xmax, Ymin, Ymax])
        - coarsetime : output one timestep every X (coarsetime=X, with X an integer)
        - singletime : output only timestep X (singletime=X, with X an integer or -1 for last)

    TODO: - make time easily accessible

    Basile de Fleurian:
    '''
    #verbosity of the code, 0 is no messages, 5 is chatty
    verbose = 0

    #first check if the user asked for some options to be applied
    for key in kwargs.keys():
        if key not in ['clipping', 'coarsetime', 'singletime']:
            raise BadOption('Provided option "{}" is not supported possibilities are : {}'.format(key, ['cliping', 'coarsetime', 'singletime']))

    if 'coarsetime' in kwargs.keys() and 'singletime' in kwargs.keys():
        raise BadOption("You can't specify both 'coarsetime' and 'singletime'")

    # File checking and creation {{{
    Dir = path.basename(filename)
    if path.exists(filename):
        print(('File {} allready exist'.format(filename)))
        newname = input('Give a new name or "delete" to replace: ')
        if newname == 'delete':
            filelist = glob(filename + '/* ')
            for oldfile in filelist:
                remove(oldfile)
        else:
            print(('New file name is {}'.format(newname)))
            filename = newname
            mkdir(filename)
    else:
        mkdir(filename)

    # }}}

    # make an alias for results {{{
    if verbose > 3:
        print('Getting accessory variables')
    res_struct = md.results
    moving_mesh = False
    if(type(res_struct) != list):
        #Getting all the solutions of the md
        solnames = dict.keys(res_struct.__dict__)
        num_of_timesteps = 1
        #%building solutionstructure
        for solution in solnames:
            #looking for multiple time steps
            try:
                if len(res_struct.__dict__[solution]) > num_of_timesteps:
                    num_of_timesteps = len(res_struct.__dict__[solution])
                    num_of_timesteps = int(num_of_timesteps)
                    #If Suface is in the resluts we considet that we have a moving mesh
                    if 'Surface' in dict.keys(res_struct.__dict__[solution][0].__dict__):
                        moving_mesh = True
            except TypeError:
                continue
    else:
        num_of_timesteps = 1
    # }}}

    # get the mesh related variables {{{
    if verbose > 3:
        print('Now treating  the mesh')
    #first get the general things
    dim = int(md.mesh.domaintype()[0])
    every_nodes = md.mesh.numberofvertices
    every_cells = md.mesh.numberofelements
    try:
        every_edges = md.mesh.numberofedges
    except AttributeError:
        #3D meshes do not have edges
        every_edges = 0

    if np.shape(md.mesh.elements)[1] == 3 or enveloppe:
        point_per_elt = 3
        celltype = 5  #triangles
    elif np.shape(md.mesh.elements)[1] == 6:
        point_per_elt = 6
        celltype = 13  #wedges
    else:
        raise BadDimension('exportVTU does not support your element type')

    #only keep the envelope and not the bulk of the results.
    if enveloppe:  #Treating enveloppe{{{
        if dim == 3:
            mesh_alti = '0'
            is_enveloppe = np.logical_or(md.mesh.vertexonbase, md.mesh.vertexonsurface)
            enveloppe_index = np.where(is_enveloppe)[0]
            convert_index = np.nan * np.ones(np.shape(md.mesh.x))
            convert_index = np.asarray([[i, np.where(enveloppe_index == i)[0][0]] for i, val in enumerate(convert_index) if any(enveloppe_index == i)])

            num_of_points = np.size(enveloppe_index)
            points = np.column_stack((md.mesh.x[enveloppe_index],
                                      md.mesh.y[enveloppe_index],
                                      md.mesh.z[enveloppe_index]))

            num_of_elt = np.size(np.where(np.isnan(md.mesh.lowerelements))) + np.size(np.where(np.isnan(md.mesh.upperelements)))
            connect = md.mesh.elements[np.where(is_enveloppe[md.mesh.elements - 1])].reshape(int(num_of_elt), 3) - 1
            for elt in range(0, num_of_elt):
                connect[elt, 0] = convert_index[np.where(convert_index == connect[elt, 0])[0], 1][0]
                connect[elt, 1] = convert_index[np.where(convert_index == connect[elt, 1])[0], 1][0]
                connect[elt, 2] = convert_index[np.where(convert_index == connect[elt, 2])[0], 1][0]

            num_of_edges = every_edges  #looks like edges is only defined on the 2d mesh
            if num_of_edges > 0:
                edges = md.mesh.edges[:, 0:2].reshape(int(num_of_edges), 2) - 1

        else:
            raise BadDimension("exportVTU can't get an enveloppe for  dimension {}".format(dim))
    # }}}

    else:  #treating mesh{{{
        #we get all the mesh, mainly defining dummies
        num_of_elt = every_cells
        connect = md.mesh.elements - 1
        num_of_edges = every_edges
        if num_of_edges > 0:
            edges = md.mesh.edges[:, 0:2].reshape(int(num_of_edges), 2) - 1
        enveloppe_index = np.arange(0, np.size(md.mesh.x))
        num_of_points = every_nodes
        if dim == 2:
            mesh_alti = input('''This is a 2D model, what should be the 3rd dimension of the mesh :
                                        1 : md.geometry.surface
                                        2 : md.geometry.base
                                        3 : md.geometry.bed
                                        4 : 0
                                        5 : Custom\n''')
            if mesh_alti == '1':
                points = np.column_stack((md.mesh.x, md.mesh.y, md.geometry.surface))
            elif mesh_alti == '2':
                points = np.column_stack((md.mesh.x, md.mesh.y, md.geometry.base))
            elif mesh_alti == '3':
                points = np.column_stack((md.mesh.x, md.mesh.y, md.geometry.bed))
            elif mesh_alti == '4':
                points = np.column_stack((md.mesh.x, md.mesh.y, 0. * md.mesh.x))
            elif mesh_alti == '5':
                alti_field = input("Which field should be used as 3rd dimension: ")
                alti_var = eval(alti_field)
                if np.shape(np.squeeze(alti_var)) == np.shape(md.mesh.x):
                    points = np.column_stack((md.mesh.x, md.mesh.y, np.squeeze(alti_var)))
                else:
                    raise BadDimension('field given for 3rd dimension should be defined on vertices {} is not.'.format(alti_field))
            else:
                points = np.column_stack((md.mesh.x, md.mesh.y, md.geometry.surface))
        elif dim == 3:
            mesh_alti = '0'
            points = np.column_stack((md.mesh.x, md.mesh.y, md.mesh.z))
        else:
            raise BadDimension('exportVTU does not support dimension {}'.format(dim))
    # }}}

    if 'clipping' in kwargs.keys():
        if kwargs['clipping'] is not None:
            # first get the boundaries and check them
            [Xmin, Xmax, Ymin, Ymax] = kwargs['clipping']
            if Xmin > Xmax:
                raise ClipError('Xmax ({}) should be larger than Xmin ({})'.format(Xmax, Xmin))
            if Ymin > Ymax:
                raise ClipError('Ymax ({}) should be larger than Ymin ({})'.format(Ymax, Ymin))
            if Xmin > np.nanmax(points[:, 0]) or Xmax < np.nanmin(points[:, 0]):
                raise ClipError('Your X boundaries [{}, {}] are outside of the model domain [{},{}]'.format(Xmin, Xmax, np.nanmin(points[:, 0]), np.nanmax(points[:, 0])))
            if Ymin > np.nanmax(points[:, 1]) or Ymax < np.nanmin(points[:, 1]):
                raise ClipError('Your Y boundaries [{}, {}] are outside of the model domain [{},{}]'.format(Ymin, Ymax, np.nanmin(points[:, 1]), np.nanmax(points[:, 1])))

            #boundaries should be fine lets do stuff
            InX = np.where(np.logical_and(points[:, 0] >= Xmin, points[:, 0] <= Xmax))
            InY = np.where(np.logical_and(points[:, 1] >= Ymin, points[:, 1] <= Ymax))

            Isinside = np.zeros(np.shape(points)[0], dtype=bool)
            clip_convert_index = np.nan * np.ones(np.shape(points)[0])

            #define the vertices that are within clipping window
            Inclipping = np.intersect1d(InX, InY)
            Isinside[Inclipping] = True
            points = points[Inclipping, :]
            num_of_points = np.shape(points)[0]

            #go thorough the elements and keep those for which one node is in the clipped arrea
            clipconnect = np.asarray([], dtype=int)
            for elt in connect:
                if set(elt).issubset(Inclipping):
                    clipconnect = np.append(clipconnect, elt, axis=0)

            #reshape
            num_of_elt = int(np.size(clipconnect) / 3)
            connect = clipconnect.reshape(num_of_elt, 3)

            clip_convert_index = np.asarray([[i, np.where(Inclipping == i)[0][0]] for i, val in enumerate(clip_convert_index) if any(Inclipping == i)])
            enveloppe_index = enveloppe_index[clip_convert_index[:, 0]]

            #convert indexing and exclude elements that are partly outside of the region
            for elt in range(0, num_of_elt):
                try:
                    connect[elt, 0] = clip_convert_index[np.where(clip_convert_index == connect[elt, 0])[0], 1][0]
                except IndexError:
                    connect[elt, 0] = -1
                try:
                    connect[elt, 1] = clip_convert_index[np.where(clip_convert_index == connect[elt, 1])[0], 1][0]
                except IndexError:
                    connect[elt, 1] = -1
                try:
                    connect[elt, 2] = clip_convert_index[np.where(clip_convert_index == connect[elt, 2])[0], 1][0]
                except IndexError:
                    connect[elt, 2] = -1

            connect = connect[np.where(connect != -1)[0], :]
            num_of_elt = np.shape(connect)[0]

            if num_of_edges > 0:
                clipedges = np.asarray([], dtype=int)
                for edge in edges:
                    if set(edge).issubset(Inclipping):
                        clipedges = np.append(clipedges, edge, axis=0)

                num_of_edges = int(np.size(clipedges) / 2)
                edges = clipedges.reshape(num_of_edges, 2)

                for edge in range(0, num_of_edges):
                    try:
                        edges[edge, 0] = clip_convert_index[np.where(clip_convert_index == edges[edge, 0])[0], 1][0]
                    except IndexError:
                        edges[edge, 0] = -1
                    try:
                        edges[edge, 1] = clip_convert_index[np.where(clip_convert_index == edges[edge, 1])[0], 1][0]
                    except IndexError:
                        edges[edge, 1] = -1
                edges = edges[np.where(edges != -1)[0], :]
                num_of_edges = np.shape(edges)[0]

    # }}}

    # write header and mesh {{{
    if verbose > 3:
        print('Now starting to write stuff')

    if 'coarsetime' in kwargs.keys():
        steplist = range(0, num_of_timesteps, kwargs['coarsetime'])
    elif 'singletime' in kwargs.keys():
        steplist = [kwargs['singletime']]
    else:
        steplist = range(0, num_of_timesteps)

    for step in steplist:
        if verbose > 2:
            print('Writing for step {}'.format(step))

        with open(('{}/{}_{}.vtu').format(filename, Dir, step), 'w+') as fid:
            fid.write('<?xml version="1.0"?>\n')
            fid.write('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">\n')
            fid.write('  <UnstructuredGrid>\n')
            fid.write('    <Piece NumberOfPoints="{}"  NumberOfCells="{}">\n'.format(num_of_points, num_of_elt + num_of_edges))
            tensors = []
            vectors = []
            scalars = []
            for sol in solnames:
                #getting the  fields in the solution
                if type(res_struct.__dict__[sol]).__name__ == 'solution':
                    spe_res_struct = res_struct.__dict__[sol].__getitem__(0)
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))
                elif type(res_struct.__dict__[sol]).__name__ in ['solutionstep', 'results']:
                    spe_res_struct = res_struct.__dict__[sol]
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))
                else:
                    print("WARNING, solution type '{}' is not recognise, exported results might be wrong".format(type(res_struct.__dict__[sol])))
                    spe_res_struct = res_struct.__dict__[sol]
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))

                loctensors, locvectors, locscalars = SortFields(fieldnames)
                tensors.extend(loctensors)
                vectors.extend(locvectors)
                scalars.extend(locscalars)
            for other in args:
                other_struct = md.__dict__[other]
                othernames = list(dict.keys(other_struct.__dict__))

                loctensors, locvectors, locscalars = SortFields(othernames)
                tensors.extend(loctensors)
                vectors.extend(locvectors)
                scalars.extend(locscalars)

            fid.write('      <PointData Scalars="{}"'.format(scalars))
            if len(vectors) > 0:
                fid.write(' Vectors="{}"'.format(vectors[:-1]))
            if len(tensors) > 0:
                fid.write(' Tensors="{}"'.format(tensors[:-2]))
            fid.write('>\n')

            saved_cells = {}
            saved_edges = {}
            saved_const = {}
            timestep = step

            # }}}
            # {{{loop over the different solution structures
            # first check if there are solutions to grab
            for sol in solnames:
                treated_res = []
                #dealing with results on different timesteps
                try:
                    if(len(res_struct.__dict__[sol]) > timestep):
                        timestep = step
                    else:
                        timestep = np.size(res_struct.__dict__[sol])
                except TypeError:
                    #result as no len() so no timesteps
                    timestep = 1

                #getting the  fields in the solution
                if(type(res_struct.__dict__[sol]).__name__ == 'solution'):
                    spe_res_struct = res_struct.__dict__[sol].__getitem__(timestep)
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))
                elif(type(res_struct.__dict__[sol]).__name__ == 'solutionstep'):
                    spe_res_struct = res_struct.__dict__[sol]
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))
                elif(type(res_struct.__dict__[sol]).__name__ == 'results'):  #this is a result without steps
                    spe_res_struct = res_struct.__dict__[sol]
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))
                else:
                    print("WARNING, solution type '{}' is not recognise, exported results might be wrong".format(type(res_struct.__dict__[sol])))
                    spe_res_struct = res_struct.__dict__[sol]
                    fieldnames = list(dict.keys(spe_res_struct.__dict__))

                tensors, vectors, ScalarNames = SortFields(fieldnames)

                #check which field is a real result and print
                for field in fieldnames:
                    if field in treated_res:
                        if verbose > 2:
                            print("{}.{} is already done".format(sol, field))
                        continue

                    elif field in vectors:
                        if verbose > 2:
                            print("Treating {}.{} as a vector ".format(sol, field))
                        TreatVector(fid, fmtout, spe_res_struct, sol, field, treated_res, enveloppe_index)

                    elif field in tensors:
                        if verbose > 2:
                            print("Treating {}.{} as a tensor ".format(sol, field))
                        TreatTensor(fid, fmtout, spe_res_struct, sol, field, treated_res, enveloppe_index)

                    else:
                        if np.size(spe_res_struct.__dict__[field]) == 1:
                            if verbose > 2:
                                print("Treating {}.{} as a constant ".format(sol, field))
                            if field == 'time':
                                current_time = spe_res_struct.__dict__[field]
                            saved_const[".".join((sol, field))] = np.squeeze(spe_res_struct.__dict__[field])

                        elif np.size(spe_res_struct.__dict__[field]) == every_nodes:
                            if verbose > 2:
                                print("Treating {}.{} as a node variable ".format(sol, field))
                            TreatScalar(fid, fmtout, spe_res_struct, sol, field, enveloppe_index)

                        elif np.shape(spe_res_struct.__dict__[field])[0] == np.size(spe_res_struct.__dict__[field]) == every_cells:
                            saved_cells[".".join((sol, field))] = np.squeeze(spe_res_struct.__dict__[field])

                        elif np.shape(spe_res_struct.__dict__[field])[0] == np.size(spe_res_struct.__dict__[field]) == every_edges and num_of_edges > 0:
                            saved_edges[".".join((sol, field))] = np.squeeze(spe_res_struct.__dict__[field])

                        else:
                            print("format for field {}.{} is not suported, field is skipped".format(sol, field))
            # }}}
            # loop on arguments, if something other than result is asked, do it now {{{
            for other in args:
                treated_res = []
                if verbose > 3:
                    print("Now treating {}".format(other))
                other_struct = md.__dict__[other]
                othernames = list(dict.keys(other_struct.__dict__))
                tensors, vectors, ScalarNames = SortFields(othernames)
                for field in othernames:
                    if field in treated_res:
                        if verbose > 2:
                            print("{}.{} is already done".format(other, field))
                        continue
                    elif field in vectors:
                        TreatVector(fid, fmtout, other_struct, other, field, treated_res, enveloppe_index)

                    elif field in tensors:
                        if verbose > 2:
                            print("Treating {}.{} as a tensor ".format(sol, field))
                        TreatTensor(fid, fmtout, other_struct, other, field, treated_res, enveloppe_index)
                        #now treating fields that are not vectors or tensors

                    else:
                        if np.size(other_struct.__dict__[field]) == 1:
                            if verbose > 2:
                                print("Treating {}.{} as an constant ".format(other, field))
                            if field == 'time':
                                current_time = other_struct.__dict__[field]
                            saved_const[".".join((other, field))] = np.squeeze(other_struct.__dict__[field])

                        elif np.size(other_struct.__dict__[field]) == every_nodes:
                            if verbose > 2:
                                print("Treating {}.{} as a node variable ".format(other, field))
                            TreatScalar(fid, fmtout, other_struct, other, field, enveloppe_index)

                        elif np.shape(other_struct.__dict__[field])[0] == every_nodes + 1:
                            if verbose > 3:
                                print("Treating {}.{} as a node forcing variable".format(other, field))
                            TreatForcing(fid, fmtout, other_struct, other, field, treated_res, enveloppe_index, current_time)

                        elif np.shape(other_struct.__dict__[field])[0] == np.size(other_struct.__dict__[field]) == every_cells:
                            if verbose > 3:
                                print("Treating {}.{} as a cell variable".format(other, field))
                            saved_cells[".".join((other, field))] = np.squeeze(other_struct.__dict__[field])

                        elif np.shape(other_struct.__dict__[field])[0] == np.size(other_struct.__dict__[field]) == every_edges and num_of_edges > 0:
                            if verbose > 3:
                                print("Treating {}.{} as an edge variable".format(other, field))
                            saved_edges[".".join((other, field))] = np.squeeze(other_struct.__dict__[field])

                        else:
                            print("format for field {}.{} is not suported, field is skipped".format(other, field))
            fid.write('      </PointData>\n')
            # }}}
            # Now writting cell variables {{{
            if np.size(list(saved_cells.keys())) > 0 or np.size(list(saved_edges.keys())) > 0:
                cellkeys = list(saved_cells.keys())
                edgekeys = list(saved_edges.keys())
                if len(cellkeys) > 0 and len(edgekeys) > 0:
                    savekeys = list(saved_cells.keys())
                    savekeys.extend(edgekeys)
                elif len(cellkeys) > 0:
                    savekeys = cellkeys
                elif len(edgekeys) > 0:
                    savekeys = edgekeys
                if verbose > 3:
                    print("Saving cell for {}".format(savekeys))
                fid.write('      <CellData Scalars="{}">\n'.format(savekeys))

            if np.size(list(saved_cells.keys())) > 0:
                for key in cellkeys:
                    outval = saved_cells[key]
                    if num_of_edges > 0:
                        if fmtout == "binary":
                            outval = np.append(outval, np.nan * np.ones((num_of_edges)))
                        else:
                            outval = np.append(outval, -9999.999 * np.ones((num_of_edges)))
                    if verbose > 3:
                        print("writing {} values of type {} for {}".format(len(outval), outval.dtype, key))

                    fid.write('        <DataArray type="Float32" Name="{}" format="{}">\n'.format(key, fmtout))
                    WriteIt(outval, fid, fmtout)
                    fid.write('        </DataArray>\n')

            # }}}
            # Now writting edge variables {{{
            if np.size(list(saved_edges.keys())) > 0:
                for key in list(saved_edges.keys()):
                    if fmtout == "binary":
                        outval = np.nan * np.ones((num_of_elt))
                    else:
                        outval = -9999.999 * np.ones((num_of_elt))
                    outval = np.append(outval, saved_edges[key])
                    fid.write('        <DataArray type="Float32" Name="{}" format="{}">\n'.format(key, fmtout))
                    WriteIt(outval, fid, fmtout)
                    fid.write('        </DataArray>\n')
            if np.size(list(saved_cells.keys())) > 0 or np.size(list(saved_edges.keys())) > 0:
                fid.write('      </CellData>\n')
            # }}}

            # Now writting constants # {{{
            if np.size(list(saved_const.keys())) > 0:
                fid.write('      <FieldData>\n')
                for key in list(saved_const.keys()):
                    fid.write('        <DataArray type="Float32" Name="{}" format="{}">\n'.format(key, fmtout))
                    WriteIt(saved_const[key], fid, fmtout)
                    fid.write('        </DataArray>\n')
                fid.write('      </FieldData>\n')
            # }}}

            #Mesh Treatment and write, it needs to loop to allow variable geometry {{{
            #updating z for mesh evolution
            if moving_mesh and mesh_alti == '1':
                points[:, 2] = np.squeeze(res_struct.__dict__['TransientSolution'][step].__dict__['Surface'][enveloppe_index])
            elif moving_mesh and mesh_alti == '2':
                points[:, 2] = np.squeeze(res_struct.__dict__['TransientSolution'][step].__dict__['Base'][enveloppe_index])

            #Now write points locations
            fid.write('      <Points>\n')
            fid.write('        <DataArray type="Float32" Name="Points" NumberOfComponents="3" format="{}">\n'.format(fmtout))
            WriteIt(points, fid, fmtout)
            fid.write('        </DataArray>\n')
            fid.write('      </Points>\n')

            #cells are a combination of element and edges
            # we need node conectivity offsets and types
            #offsets is the cummulative index of the last elemant of each cell (1 indexed)
            flat_elt = connect.flatten()
            elt_offset = np.arange(0, num_of_elt * point_per_elt, point_per_elt, dtype=np.int64) + point_per_elt
            elt_type = celltype * np.ones((num_of_elt), dtype=np.uint8)
            if num_of_edges > 0:
                flat_edges = edges.flatten()
                flat_cells = np.hstack((flat_elt, flat_edges))
                edge_offset = np.arange(0, num_of_edges * 2, 2) + 2 + elt_offset[-1]
                cell_offset = np.hstack((elt_offset, edge_offset))
                edge_type = 3 * np.ones((num_of_edges), dtype=np.uint8)
                cell_type = np.hstack((elt_type, edge_type))
            else:
                flat_cells = flat_elt
                cell_offset = elt_offset
                cell_type = elt_type

            if verbose > 3:
                print("""writing mesh structure:
                                  connectivity of shape {}
                                  cell offset of shape {}
                                  cell types of shape{}""".format(np.shape(flat_cells), np.shape(cell_offset), np.shape(cell_type)))
            #write cells Informations
            fid.write('      <Cells>\n')
            fid.write('        <DataArray type="Int64" Name="connectivity" format="{}">\n'.format(fmtout))
            WriteIt(flat_cells, fid, fmtout)
            fid.write('        </DataArray>\n')
            fid.write('        <DataArray type="Int64" Name="offsets" format="{}">\n'.format(fmtout))
            WriteIt(cell_offset, fid, fmtout)
            fid.write('        </DataArray>\n')
            fid.write('        <DataArray type="UInt8" Name="types" format="{}">\n'.format(fmtout))
            WriteIt(cell_type, fid, fmtout)
            fid.write('        </DataArray>\n')
            fid.write('      </Cells>\n')
            fid.write('    </Piece>\n')
            fid.write('  </UnstructuredGrid>\n')
            fid.write('</VTKFile>\n')
            # }}}


def SortFields(fieldnames):
    #we check on sizes so there is a slight chance that logs can be picked as results, we remove them to avoid that
    for trashfield in ['errlog', 'outlog']:
        if trashfield in fieldnames:
            fieldnames.remove(trashfield)

    #Sorting scalars, vectors and tensors
    tensors = [field for field in fieldnames if field[-2:] in ['xx', 'yy', 'xy', 'zz', 'xz', 'yz']]
    non_tensor = [field for field in fieldnames if field not in tensors]
    vectors = [field for field in non_tensor if field[-1] in ['x', 'y', 'z']]
    #get the name of scalar fields remove, vectors, tensors and things that are not proper results
    scalars = [field for field in fieldnames if field not in tensors + vectors]
    dump = ["ConvergenceNumSteps", "step", "time"]
    for trash in dump:
        try:
            scalars.remove(trash)
        except ValueError:
            [scalars.remove(name) for name in scalars if trash in name]
            continue
    #clean up vector and tensors that might be here and should not
    # we check that at least two of the vector component are here
    for namelist in [vectors, tensors]:
        for name in list(namelist):
            coord = name[-1]
            if coord == 'x' and name[:-1] + 'y' in namelist:
                continue
            elif coord == 'y' and name[:-1] + 'x' in namelist:
                continue
            elif coord == 'z' and name[:-1] + 'x' in namelist:
                continue
            else:
                scalars.extend([name])
                namelist.remove(name)
    return tensors, vectors, scalars


def TreatScalar(fid, fmtout, structure, structname, fieldname, enveloppe_index):
    array = np.squeeze(structure.__dict__[fieldname][enveloppe_index])
    fid.write('        <DataArray type="Float32" Name="{}" NumberOfComponents="1" format="{}">\n'.format(".".join((structname, fieldname)), fmtout))
    WriteIt(array, fid, fmtout)
    fid.write('        </DataArray>\n')


def TreatVector(fid, fmtout, structure, structname, fieldname, treated_res, enveloppe_index):
    Vxstruct = np.squeeze(structure.__dict__[fieldname[:-1] + 'x'])
    Vystruct = np.squeeze(structure.__dict__[fieldname[:-1] + 'y'])
    Vx = Vxstruct[enveloppe_index]
    Vy = Vystruct[enveloppe_index]
    treated_res += [fieldname[:-1] + 'x', fieldname[:-1] + 'y']
    try:
        Vzstruct = np.squeeze(structure.__dict__[fieldname[:-1] + 'z'])
        treated_res += [fieldname[:-1] + 'z']
        Vz = Vzstruct[enveloppe_index]
    except KeyError:
        Vz = np.zeros(np.shape(Vx))
    Vector = (np.vstack((Vx, Vy, Vz)).T).flatten()
    fid.write('        <DataArray type="Float32" Name="{}" NumberOfComponents="3" format="{}">\n'.format(".".join((structname, fieldname[:-1])), fmtout))
    WriteIt(Vector, fid, fmtout)
    fid.write('        </DataArray>\n')


def TreatTensor(fid, fmtout, structure, structname, fieldname, treated_res, enveloppe_index):
    Txxstruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'xx'])
    Txystruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'xy'])
    Tyystruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'yy'])
    treated_res += [fieldname[:-2] + 'xx', fieldname[:-2] + 'xy', fieldname[:-2] + 'yy']
    Txx = Txxstruct[enveloppe_index]
    Tyy = Tyystruct[enveloppe_index]
    Txy = Txystruct[enveloppe_index]
    try:
        Tzzstruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'zz'])
        Txzstruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'xz'])
        Tyzstruct = np.squeeze(structure.__dict__[fieldname[:-2] + 'yz'])
        treated_res += [fieldname[:-2] + 'zz', fieldname[:-2] + 'xz', fieldname[:-2] + 'yz']
        Tzz = Tzzstruct[enveloppe_index]
        Txz = Txzstruct[enveloppe_index]
        Tyz = Tyzstruct[enveloppe_index]
    except KeyError:
        Tzz = np.zeros(np.shape(Txx))
        Txz = np.zeros(np.shape(Txx))
        Tyz = np.zeros(np.shape(Txx))

    Tensor = (np.vstack((Txx, Tyy, Tzz, Txy, Tyz, Txz)).T).flatten()
    fid.write('        <DataArray type="Float32" Name="{}" NumberOfComponents="6" format="{}">\n'.format(".".join((structname, fieldname[:-1])), fmtout))
    WriteIt(Tensor, fid, fmtout)
    fid.write('        </DataArray>\n')


def TreatForcing(fid, fmtout, structure, structname, fieldname, treated_res, enveloppe_index, current_time):
    #we are dealing with a forcing of some kind.
    forcing_time = structure.__dict__[fieldname][-1, :]
    if any(forcing_time == current_time):
        forcing_index = np.where(forcing_time == current_time)
        forcing_val = structure.__dict__[fieldname][:, forcing_index]
    elif forcing_time[0] > current_time:
        forcing_val = structure.__dict__[fieldname][:, 0]
    elif forcing_time[-1] < current_time:
        forcing_val = structure.__dict__[fieldname][:, -1]
    else:
        forcing_index = np.where(forcing_time < current_time)[-1][-1]
        delta_time = forcing_time[forcing_index + 1] - forcing_time[forcing_index]  #compute forcing Dt
        delta_current = current_time - forcing_time[forcing_index]  # time since last forcing
        ratio = delta_current / delta_time  #compute weighting factor for preceding forcing vallue
        forcing_evol = (structure.__dict__[fieldname][:, forcing_index + 1] - structure.__dict__[fieldname][:, forcing_index]) * ratio
        forcing_val = structure.__dict__[fieldname][:, forcing_index] + forcing_evol
    array = forcing_val[enveloppe_index]
    # and now write it down
    fid.write('        <DataArray type="Float32" Name="{}" NumberOfComponents="1" format="{}">\n'.format(".".join((structname, fieldname)), fmtout))
    WriteIt(array, fid, fmtout)
    fid.write('        </DataArray>\n')


def WriteIt(Data, fid, fmtout):
    vtu_to_numpy_type = {
        "Float32": np.dtype(np.float32),
        "Float64": np.dtype(np.float64),
        "Int8": np.dtype(np.int8),
        "Int16": np.dtype(np.int16),
        "Int32": np.dtype(np.int32),
        "Int64": np.dtype(np.int64),
        "UInt8": np.dtype(np.uint8),
        "UInt16": np.dtype(np.uint16),
        "UInt32": np.dtype(np.uint32),
        "UInt64": np.dtype(np.uint64),
    }
    if fmtout == 'binary':
        try:
            datatype = Data.dtype
        except AttributeError:
            datatype = type(Data)
        if datatype == np.float64:
            Data = np.float32(Data)
        try:
            data_bytes = Data.tobytes()
        except AttributeError:
            data_bytes = np.asarray(Data).tobytes()
        # collect header
        header = np.array(len(data_bytes), dtype=vtu_to_numpy_type['UInt32'])
        fid.write(b64encode(header.tobytes() + data_bytes).decode())
        fid.write('\n')
        #cell_type.tofile(fid)
    elif fmtout == 'ascii':
        np.savetxt(fid, Data, fmt='%g')


def cleanOutliers(Val, fmtout):
    #paraview does not like NaN in ascii files, replacing
    if np.isnan(Val):
        if fmtout == 'ascii':
            CleanVal = -9999.999

    #also checking for very small value that mess up
    elif (abs(Val) < 1.0e-20):
        CleanVal = 0.0
    else:
        CleanVal = Val
    return CleanVal


class BadDimension(Exception):
    """The required dimension is not supported yet."""


class BadOption(Exception):
    """The given option does not exist."""


class ClipError(Exception):
    """Error while trying to clip the domain."""
