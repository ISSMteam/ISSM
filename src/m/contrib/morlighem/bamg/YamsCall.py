import numpy as np
import MatlabFuncs as m
import time
import subprocess
import os
from ComputeHessian import ComputeHessian
from ComputeMetric import ComputeMetric


def YamsCall(md, field, hmin, hmax, gradation, epsilon):
    """
    YAMSCALL - call yams

       build a metric using the Hessian of the given field
       call Yams and the output mesh is plugged onto the model
     - hmin = minimum edge length (m)
     - hmax = maximum edge length (m)
     - gradation = maximum edge length gradation between 2 elements
     - epsilon = average error on each element (m / yr)

       Usage:
          md = YamsCall(md, field, hmin, hmax, gradation, epsilon)

       Example:
          md = YamsCall(md, md.inversion.vel_obs, 1500, 1.0e8, 1.3, 0.9)
    """

    #2d geometric parameter (do not change)
    scale = 2. / 9.

    #Compute Hessian
    t1 = time.time()
    print(("%s" % '      computing Hessian...'))
    hessian = ComputeHessian(md.mesh.elements, md.mesh.x, md.mesh.y, field, 'node')
    t2 = time.time()
    print(("%s%d%s\n" % (' done (', t2 - t1, ' seconds)')))

    #Compute metric
    t1 = time.time()
    print(("%s" % '      computing metric...'))
    metric = ComputeMetric(hessian, scale, epsilon, hmin, hmax, np.empty(0, int))
    t2 = time.time()
    print(("%s%d%s\n" % (' done (', t2 - t1, ' seconds)')))

    #write files
    t1 = time.time()
    print(("%s" % '      writing initial mesh files...'))
    np.savetxt('carre0.met', metric)

    f = open('carre0.mesh', 'w')

    #initialiation
    f.write("\n%s\n%i\n" % ('MeshVersionFormatted', 1))

    #dimension
    f.write("\n%s\n%i\n" % ('Dimension', 2))

    #Vertices
    f.write("\n%s\n%i\n\n" % ('Vertices', md.mesh.numberofvertices))
    for i in range(0, md.mesh.numberofvertices):
        f.write("%8g %8g %i\n" % (md.mesh.x[i], md.mesh.y[i], 0))

    #Triangles
    f.write("\n\n%s\n%i\n\n" % ('Triangles', md.mesh.numberofelements))
    for i in range(0, md.mesh.numberofelements):
        f.write("%i %i %i %i\n" % (md.mesh.elements[i, 0], md.mesh.elements[i, 1], md.mesh.elements[i, 2], 0))
    numberofelements1 = md.mesh.numberofelements

    #Deal with rifts
    if np.any(not np.isnan(md.rifts.riftstruct)):

        #we have the list of triangles that make up the rift. keep those triangles around during refinement.
        triangles = np.empty(0, int)
        for riftstruct in md.rifts.riftstruct:
            triangles = np.concatenate((triangles, riftstruct.segments[:, 2]))

        f.write("\n\n%s\n%i\n\n" % ('RequiredTriangles', np.size(triangles)))
        for triangle in triangles:
            f.write("%i\n" % triangle)

    #close
    f.close()
    t2 = time.time()
    print(("%s%d%s\n" % (' done (', t2 - t1, ' seconds)')))

    #call yams
    print(("%s\n" % '      call Yams...'))
    if m.ispc():
        #windows
        subprocess.call('yams2 - win - O 1 - v - 0 - ecp - hgrad %g carre0 carre1' % gradation, shell=True)
    elif m.ismac():
        #Macosx
        subprocess.call('yams2 - osx - O 1 - v - 0 - ecp - hgrad %g carre0 carre1' % gradation, shell=True)
    else:
        #Linux
        subprocess.call('yams2 - linux - O 1 - v - 0 - ecp - hgrad %g carre0 carre1' % gradation, shell=True)

    #plug new mesh
    t1 = time.time()
    print(("\n%s" % '      reading final mesh files...'))
    Tria = np.loadtxt('carre1.tria', int)
    Coor = np.loadtxt('carre1.coor', float)
    md.mesh.x = Coor[:, 0]
    md.mesh.y = Coor[:, 1]
    md.mesh.z = np.zeros((np.size(Coor, axis=0), 1))
    md.mesh.elements = Tria
    md.mesh.numberofvertices = np.size(Coor, axis=0)
    md.mesh.numberofelements = np.size(Tria, axis=0)
    numberofelements2 = md.mesh.numberofelements
    t2 = time.time()
    print(("%s%d%s\n\n" % (' done (', t2 - t1, ' seconds)')))

    #display number of elements
    print(("\n%s %i" % ('      inital number of elements:', numberofelements1)))
    print(("\n%s %i\n\n" % ('      new    number of elements:', numberofelements2)))

    #clean up:
    os.remove('carre0.mesh')
    os.remove('carre0.met')
    os.remove('carre1.tria')
    os.remove('carre1.coor')
    os.remove('carre1.meshb')

    return md
