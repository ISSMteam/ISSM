import numpy as np
from os import path, remove, mkdir
from glob import glob


def exportVTK(filename, md, *args, enveloppe=False, **kwargs):
    '''
    vtk export
    function exportVTK(filename, md)
    creates a directory with the vtk files for displays in paraview
    (only work for triangle and wedges based on their number of nodes)

    Usage:
    exportVTK('DirName', md)
    exportVTK('DirName', md, 'geometry', 'mesh')
    exportVTK('DirName', md, 'geometry', 'mesh', enveloppe = True)

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

    print("""
    =========================================
    #     A                                 #
    #    / \      exportVTK is now obsolete #
    #   / | \     You should use export VTU #
    #  /  |  \    faster, smaller files     #
    # /   o   \   and more capacities       #
    # ---------                             #
    #========================================
    """)


    for key in kwargs.keys():
        if key not in ['clipping', 'coarsetime', 'singletime']:
            raise BadOption('Provided option "{}" is not supported possibilities are : {}'.format(key, ['cliping', 'coarsetime', 'singletime']))

    if 'coarsetime' in kwargs.keys() and 'singletime' in kwargs.keys():
        raise BadOption("You can't specify both 'coarsetime' and 'singletime'")

    # File checking and creation {{{
    Dir = path.basename(filename)
    Path = filename[:-len(Dir)]
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

    # this is the result structure {{{
    if verbose > 3:
        print('Getting accessorie variables')
    res_struct = md.results
    moving_mesh = False
    if(type(res_struct) != list):
        #Getting all the solutions of the md
        solnames = dict.keys(res_struct.__dict__)
        num_of_sols = len(solnames)
        num_of_timesteps = 1
        #%building solutionstructure
        for solution in solnames:
            #looking for multiple time steps
            try:
                if len(res_struct.__dict__[solution]) > num_of_timesteps:
                    num_of_timesteps = len(res_struct.__dict__[solution])
                    num_of_timesteps = int(num_of_timesteps)
                    if 'Surface' in dict.keys(res_struct.__dict__[solution][0].__dict__):
                        moving_mesh = True
            except TypeError:
                continue
    else:
        num_of_timesteps = 1
    # }}}

    # get the element related variables {{{
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
        raise BadDimension('exportVTK does not support your element type')

    #only keep the envelope and not the bulk of the results.
    if enveloppe:
        if dim == 3:
            mesh_alti = '1'
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
            raise BadDimension("exportVTK can't get an enveloppe for  dimension {}".format(dim))

    else:
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
            mesh_alti = '1'
            points = np.column_stack((md.mesh.x, md.mesh.y, md.mesh.z))
        else:
            raise BadDimension('exportVTK does not support dimension {}'.format(dim))

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
        saved_cells = {}
        saved_edges = {}
        timestep = step
        with open((filename + '/Timestep.vtk' + str(timestep) + '.vtk'), 'w+') as fid:
            fid.write('# vtk DataFile Version 3.0 \n')
            fid.write('Data for run {} \n'.format(md.miscellaneous.name))
            fid.write('ASCII \n')
            fid.write('DATASET UNSTRUCTURED_GRID \n')
            fid.write('POINTS {:d} float\n'.format(num_of_points))
            #updating z for mesh evolution
            if moving_mesh and mesh_alti in ['1', '2']:
                base = np.squeeze(res_struct.__dict__['TransientSolution'][step].__dict__['Base'][enveloppe_index])
                thick_change_ratio = (np.squeeze(res_struct.__dict__['TransientSolution'][step].__dict__['Thickness'][enveloppe_index]) / md.geometry.thickness[enveloppe_index])
                above_bed = points[:, 2] - md.geometry.base[enveloppe_index]
                altitude = base + thick_change_ratio * above_bed
            else:
                altitude = points[:, 2]
            for index, point in enumerate(points):
                fid.write('{:f} {:f} {:f} \n'.format(point[0], point[1], altitude[index]))

            fid.write('CELLS {:d} {:d}\n'.format((num_of_elt + num_of_edges), num_of_elt  * (point_per_elt + 1) + num_of_edges * 3))

            for elt in range(0, num_of_elt):
                if celltype == 5:
                    fid.write('3 {:d} {:d} {:d}\n'.format(connect[elt, 0],
                                                          connect[elt, 1],
                                                          connect[elt, 2]))
                elif celltype == 13:
                    fid.write('6 {:d} {:d} {:d} {:d} {:d} {:d}\n'.format(connect[elt, 0],
                                                                         connect[elt, 1],
                                                                         connect[elt, 2],
                                                                         connect[elt, 3],
                                                                         connect[elt, 4],
                                                                         connect[elt, 5]))
            for edge in range(0, num_of_edges):
                fid.write('2 {:d} {:d}\n'.format(edges[edge, 0],
                                                 edges[edge, 1]))

            fid.write('CELL_TYPES {:d}\n'.format(num_of_elt + num_of_edges))
            for elt in range(0, num_of_elt):
                fid.write('{:d}\n'.format(celltype))
                for edge in range(0, num_of_edges):
                    fid.write('3\n')  #3 is for lines

            fid.write('POINT_DATA {:s} \n'.format(str(num_of_points)))
            # }}}
            # {{{loop over the different solution structures
            # first check if there are solutions to grab
            if 'solnames' in locals():
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

                    #Sorting scalars, vectors and tensors
                    tensors = [field for field in fieldnames if field[-2:] in ['xx', 'yy', 'xy', 'zz', 'xz', 'yz']]
                    non_tensor = [field for field in fieldnames if field not in tensors]
                    vectors = [field for field in non_tensor if field[-1] in ['x', 'y', 'z'] and field[-4:] not in ['Flux']]
                    #check which field is a real result and print
                    for field in fieldnames:
                        if verbose > 2:
                            print("Treating {}".format(field))
                        if field in treated_res:
                            if verbose > 2:
                                print("{} is already done".format(field))
                            continue
                        elif field in vectors:
                            if verbose > 2:
                                print("{} is a vector".format(field))
                            try:
                                Vxstruct = np.squeeze(spe_res_struct.__dict__[field[:-1] + 'x'])
                                Vystruct = np.squeeze(spe_res_struct.__dict__[field[:-1] + 'y'])
                                treated_res += [field[:-1] + 'x', field[:-1] + 'y']
                                if dim == 3 and field[:-1] + 'z' in fieldnames:
                                    #some fields like adjoint or always 2D
                                    Vzstruct = np.squeeze(spe_res_struct.__dict__[field[:-1] + 'z'])
                                    treated_res += [field[:-1] + 'z']

                            except KeyError:
                                fieldnames += field
                                vectors.remove(field)

                            fid.write('VECTORS {} float \n'.format(field[:-1]))
                            for node in range(0, num_of_points):
                                Vx = cleanOutliers(Vxstruct[enveloppe_index[node]])
                                Vy = cleanOutliers(Vystruct[enveloppe_index[node]])
                                if dim == 3 and field[:-1] + 'z' in fieldnames:
                                    Vz = cleanOutliers(Vzstruct[enveloppe_index[node]])
                                    fid.write('{:f} {:f} {:f}\n'.format(Vx, Vy, Vz))
                                else:
                                    fid.write('{:f} {:f} {:f}\n'.format(Vx, Vy, 0))

                        elif field in tensors:
                            if verbose > 2:
                                print("{} is a tensor".format(field))
                            try:
                                Txxstruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'xx'])
                                Txystruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'xy'])
                                Tyystruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'yy'])
                                treated_res += [field[:-2] + 'xx', field[:-2] + 'xy', field[:-2] + 'yy']
                                if dim == 3:
                                    Tzzstruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'zz'])
                                    Txzstruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'xz'])
                                    Tyzstruct = np.squeeze(spe_res_struct.__dict__[field[:-2] + 'yz'])
                                    treated_res += [field[:-2] + 'zz', field[:-2] + 'xz', field[:-2] + 'yz']

                            except KeyError:
                                fieldnames += field
                                tensors.remove(field)

                            fid.write('TENSORS {} float \n'.format(field[:-2]))
                            for node in range(0, num_of_points):
                                Txx = cleanOutliers(Txxstruct[enveloppe_index[node]])
                                Tyy = cleanOutliers(Tyystruct[enveloppe_index[node]])
                                Txy = cleanOutliers(Txystruct[enveloppe_index[node]])
                                if dim == 3:
                                    Tzz = cleanOutliers(Tzzstruct[enveloppe_index[node]])
                                    Txz = cleanOutliers(Txzstruct[enveloppe_index[node]])
                                    Tyz = cleanOutliers(Tyzstruct[enveloppe_index[node]])
                                    fid.write('{:f} {:f} {:f}\n'.format(Txx, Txy, Txz))
                                    fid.write('{:f} {:f} {:f}\n'.format(Txy, Tyy, Tyz))
                                    fid.write('{:f} {:f} {:f}\n'.format(Txz, Tyz, Tzz))
                                elif dim == 2:
                                    fid.write('{:f} {:f} {:f}\n'.format(Txx, Txy, 0))
                                    fid.write('{:f} {:f} {:f}\n'.format(Txy, Tyy, 0))
                                    fid.write('{:f} {:f} {:f}\n'.format(0, 0, 0))
                        else:
                            if np.size(spe_res_struct.__dict__[field]) == 1:
                                if field == 'time':
                                    current_time = spe_res_struct.__dict__[field]
                                    #skipping integers
                                continue
                            elif np.size(spe_res_struct.__dict__[field]) == every_nodes:
                                fid.write('SCALARS {} float 1 \n'.format(field))
                                fid.write('LOOKUP_TABLE default\n')
                                for node in range(0, num_of_points):
                                    outval = cleanOutliers(np.squeeze(spe_res_struct.__dict__[field][enveloppe_index[node]]))
                                    fid.write('{:f}\n'.format(outval))
                            elif np.shape(spe_res_struct.__dict__[field])[0] == np.size(spe_res_struct.__dict__[field]) == every_cells:
                                saved_cells[field] = np.squeeze(spe_res_struct.__dict__[field])
                            elif np.shape(spe_res_struct.__dict__[field])[0] == np.size(spe_res_struct.__dict__[field]) == every_edges:
                                saved_edges[field] = np.squeeze(spe_res_struct.__dict__[field])
                            else:
                                print("format for field {}.{} is not suported, field is skipped".format(sol, field))
            # }}}
            # loop on arguments, if something other than result is asked, do it now {{{
            for other in args:
                other_struct = md.__dict__[other]
                othernames = (dict.keys(other_struct.__dict__))
                for field in othernames:
                    if np.size(other_struct.__dict__[field]) == 1:
                        #skipping integers
                        continue
                    elif np.size(other_struct.__dict__[field]) == every_nodes:
                        fid.write('SCALARS {} float 1 \n'.format(field))
                        fid.write('LOOKUP_TABLE default\n')
                        for node in range(0, num_of_points):
                            outval = cleanOutliers(other_struct.__dict__[field][enveloppe_index[node]])
                            fid.write('{:f}\n'.format(outval))
                    elif np.shape(other_struct.__dict__[field])[0] == every_nodes + 1:
                        #we are dealing with a forcing of some kind.
                        forcing_time = other_struct.__dict__[field][-1, :]
                        if any(forcing_time == current_time):
                            forcing_index = np.where(forcing_time == current_time)
                            forcing_val = other_struct.__dict__[field][:, forcing_index]
                        elif forcing_time[0] > current_time:
                            forcing_val = other_struct.__dict__[field][:, 0]
                        elif forcing_time[-1] < current_time:
                            forcing_val = other_struct.__dict__[field][:, -1]
                        else:
                            forcing_index = np.where(forcing_time < current_time)[-1][-1]
                            delta_time = forcing_time[forcing_index + 1] - forcing_time[forcing_index]  #compute forcing Dt
                            delta_current = current_time - forcing_time[forcing_index]  # time since last forcing
                            ratio = delta_current / delta_time  #compute weighting factor for preceding forcing vallue
                            forcing_evol = (other_struct.__dict__[field][:, forcing_index + 1] - other_struct.__dict__[field][:, forcing_index]) * ratio
                            forcing_val = other_struct.__dict__[field][:, forcing_index] + forcing_evol
                        # and now write it down
                        fid.write('SCALARS {}_{} float 1 \n'.format(other, field))
                        fid.write('LOOKUP_TABLE default\n')
                        for node in range(0, num_of_points):
                            outval = cleanOutliers(forcing_val[enveloppe_index[node]])
                            fid.write('{:f}\n'.format(outval))
                    elif np.shape(other_struct.__dict__[field])[0] == np.size(other_struct.__dict__[field]) == every_cells:
                        saved_cells[field] = other_struct.__dict__[field]
                    elif np.shape(other_struct.__dict__[field])[0] == np.size(other_struct.__dict__[field]) == every_edges:
                        saved_edges[field] = other_struct.__dict__[field]
                    else:
                        print("format for field {}.{} is not suported, field is skipped".format(other, field))
                        continue
            # }}}
            # Now writting cell variables {{{
            if np.size(list(saved_cells.keys())) > 0:
                fid.write('CELL_DATA {:d} \n'.format(num_of_elt + num_of_edges))
                for key in list(saved_cells.keys()):
                    fid.write('SCALARS {} float 1 \n'.format(key))
                    fid.write('LOOKUP_TABLE default\n')
                    for cell in range(0, num_of_elt):
                        outval = cleanOutliers(saved_cells[key][cell])
                        fid.write('{:f}\n'.format(outval))
                    for edge in range(0, num_of_edges):
                        fid.write('{:f}\n'.format(-9999.999))
            # }}}
            # Now writting edge variables {{{
            if np.size(list(saved_edges.keys())) > 0:
                for key in list(saved_edges.keys()):
                    fid.write('SCALARS {} float 1 \n'.format(key))
                    fid.write('LOOKUP_TABLE default\n')
                    for cell in range(0, num_of_elt):
                        fid.write('{:f}\n'.format(-9999.999))
                    for edge in range(0, num_of_edges):
                        outval = cleanOutliers(saved_edges[key][edge])
                        fid.write('{:f}\n'.format(outval))
    # }}}


def cleanOutliers(Val):
    #paraview does not like NaN, replacing
    if np.isnan(Val):
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
