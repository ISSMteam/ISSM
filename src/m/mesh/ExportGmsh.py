import numpy as np
"""
EXPORTGMSH - export mesh to gmsh format

   Usage:
      ExportGmsh(md,filename)
"""


def ExportGmsh(md, filename):

    print('writing gmsh mesh file')
    with open(filename, 'w+') as GmshFile:

        #initialization
        GmshFile.write('$MeshFormat \n')
        GmshFile.write('2.2 0 8 \n')
        GmshFile.write('$EndMeshFormat \n')
        GmshFile.write('$Nodes \n')
        GmshFile.write('{} \n'.format(md.mesh.numberofvertices))

        #printing point positions
        for j, eX in enumerate(md.mesh.x):
            GmshFile.write('{:g} {:14.7e} {:14.7e} 0.0 \n'.format(j + 1, eX, md.mesh.y[j]))

        GmshFile.write('$EndNodes \n')
        GmshFile.write('$Elements \n')
        GmshFile.write('{:d} \n'.format(md.mesh.numberofelements + np.shape(md.mesh.segments)[0]))

        #printing elements caracteristics for boundaries
        for ind, segment in enumerate(md.mesh.segments):
            if md.mesh.x[segment[0]] == np.nanmax(md.mesh.x) and md.mesh.x[segment[1]] == np.nanmax(md.mesh.x):
                bc_id = 1
            elif md.mesh.y[segment[0]] == np.nanmax(md.mesh.y) and md.mesh.y[segment[1]] == np.nanmax(md.mesh.y):
                bc_id = 2
            elif md.mesh.x[segment[0]] == np.nanmin(md.mesh.x) and md.mesh.x[segment[1]] == np.nanmin(md.mesh.x):
                bc_id = 3
            elif md.mesh.y[segment[0]] == np.nanmin(md.mesh.y) and md.mesh.y[segment[1]] == np.nanmin(md.mesh.y):
                bc_id = 4
            else:
                bc_id = 0

            GmshFile.write('{:g} 1 2 {:g} 1 {:g} {:g} \n'.format(ind + 1, bc_id, segment[0], segment[1]))

        #and for the body
        body_id = 1
        for elt, element in enumerate(md.mesh.elements):
            GmshFile.write('{:g} 2 2 {:g} 3 {:g} {:g} {:g} \n'.format(elt + 1, body_id, element[0], element[1], element[2]))
        GmshFile.write('$EndElements \n')
        #close
        GmshFile.close()
