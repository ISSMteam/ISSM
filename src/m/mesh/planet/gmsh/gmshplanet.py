import os
import subprocess
import numpy as np
from issmdir import issmdir
from MatlabFuncs import *
from mesh3dsurface import *
from pairoptions import *


def gmshplanet(*args):
    """gmshplanet - mesh generation for a sphere. Very specific code for Gmsh from $ISSM_DIR/src/demos/simple_geo/sphere.geo

    Available options (for more details see ISSM website http://issm.jpl.nasa.gov/):
    - radius:             radius of the planet in km
    - resolution:         resolution in km
    - refine:             provide mesh
    - refinemetric:       mesh quantity to specify resolution

    Returns 'mesh3dsurface' type mesh

    Examples:
        md.mesh = gmshplanet('radius', 6000, 'resolution', 100);
        md.mesh = gmshplanet('radius', 6000, 'resolution', 100);
    """

    # Get Gmsh version
    subproc_args = 'gmsh -info 2>&1 | command grep \'Version\' | sed -e \'s/Version[[:blank:]]*:[[:blank:]]//\' | cut -d \'.\' -f1'
    subproc = subprocess.Popen(subproc_args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = subproc.communicate()
    try:
        strErrs = errs.decode()
    except AttributeError:  # this is not a byte variable, let's assume string
        strErrs = errs
    if strErrs != '':
        # gmsh executable may not be on path; attempt to find it
        paths = [
            os.environ.get('ISSM_EXT_DIR') + '/shared/gmsh/install/bin',
            os.environ.get('ISSM_EXT_DIR') + '/static/gmsh/install/bin',
            os.environ.get('ISSM_EXT_DIR') + '/gmsh/install/bin',
            issmdir() + '/externalpackages/gmsh/install/bin',
            issmdir() + '/bin',
            '/usr/bin'
        ]
        gmshpath = ''
        for path in paths:
            if exists(path + '/gmsh'):
                gmshpath = path
                break
        if gmshpath == '':
            error('gmshplanet: gmsh executable not found!')

        os.environ['PATH'] = gmshpath ':' os.environ.get['PATH']

        # Get Gmsh version
        subproc_args = 'gmsh -info 2>&1 | command grep \'Version\' | sed -e \'s/Version[[:blank:]]*:[[:blank:]]//\' | cut -d \'.\' -f1'
        subproc = subprocess.Popen(subproc_args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        outs, errs = subproc.communicate()
        try:
            strErrs = errs.decode()
        except AttributeError:  # this is not a byte variable, let's assume string
            strErrs = errs
        if strErrs != '':
            raise Exception('gmshplanet: call to gmsh failed: {}'.format(errs))

    gmshmajorversion = int(outs)
    if gmshmajorversion not in [3, 4]:
        raise RuntimeError('gmshplanet: Gmsh major version {} not supported!'.format(gmshmajorversion))

    # Process options
    options = pairoptions(*args)
    #options = deleteduplicates(options, 1)

    # Recover parameters
    radius = options.getfieldvalue('radius') * 1000
    resolution = options.getfieldvalue('resolution') * 1000

    # Initialize mesh
    mesh = mesh3dsurface()
    # Create .geo file:  {{{
    fid = open('sphere.geo', 'w')

    # Call gmsh
    #
    # NOTE:
    # - The default format in Gmsh 3 is "msh2". Rather than conditionally
    # modifying our parsing scheme for Gmsh 4, for now, we simply set the
    # 'Mesh.MshFileVersion' option.
    # - Decreasing the value of the 'Mesh.RandomFactor' option leads to an
    # equal number of nodes and elements being produced under macOS and Linux
    # at certain resolutions using certain meshing algorithms.
    #
    fid.write('Mesh.Algorithm = 1;\n') # MeshAdapt
    fid.write('Mesh.MshFileVersion = 2;\n')
    fid.write('Mesh.RandomFactor = 1e-10;\n')
    if options.exist('refine'):
        fid.write('Mesh.Algorithm = 7;\n') # BAMG
        fid.write('Mesh.CharacteristicLengthFromPoints = 0;\n')
        if gmshmajorversion == 3:
            fid.write('Mesh.RemeshAlgorithm = 1;\n')
    fid.write('resolution = %g;\n' % resolution)
    fid.write('radius = %g;\n' % radius)
    fid.write('Point(1) = {0.0, 0.0, 0.0, resolution};\n')
    fid.write('Point(2) = {radius, 0.0, 0.0, resolution};\n')
    fid.write('Point(3) = {0, radius, 0.0, resolution};\n')
    fid.write('Circle(1) = {2, 1, 3};\n')
    fid.write('Point(4) = {-radius, 0, 0.0, resolution};\n')
    fid.write('Point(5) = {0, -radius, 0.0, resolution};\n')
    fid.write('Circle(2) = {3, 1, 4};\n')
    fid.write('Circle(3) = {4, 1, 5};\n')
    fid.write('Circle(4) = {5, 1, 2};\n')
    fid.write('Point(6) = {0, 0, -radius, resolution};\n')
    fid.write('Point(7) = {0, 0, radius, resolution};\n')
    fid.write('Circle(5) = {3, 1, 6};\n')
    fid.write('Circle(6) = {6, 1, 5};\n')
    fid.write('Circle(7) = {5, 1, 7};\n')
    fid.write('Circle(8) = {7, 1, 3};\n')
    fid.write('Circle(9) = {2, 1, 7};\n')
    fid.write('Circle(10) = {7, 1, 4};\n')
    fid.write('Circle(11) = {4, 1, 6};\n')
    fid.write('Circle(12) = {6, 1, 2};\n')

    if gmshmajorversion == 3:
        curvename = 'Line Loop'
    elif gmshmajorversion == 4:
        curvename = 'Curve Loop'

    fid.write('%s(13) = {2,8,-10};\n' % curvename)
    fid.write('Surface(14) = {13};\n')
    fid.write('%s(15) = {10,3,7};\n' % curvename)
    fid.write('Surface(16) = {15};\n')
    fid.write('%s(17) = {-8,-9,1};\n' % curvename)
    fid.write('Surface(18) = {17};\n')
    fid.write('%s(19) = {-11,-2,5};\n' % curvename)
    fid.write('Surface(20) = {19};\n')
    fid.write('%s(21) = {-5,-12,-1};\n' % curvename)
    fid.write('Surface(22) = {21};\n')
    fid.write('%s(23) = {-3,11,6};\n' % curvename)
    fid.write('Surface(24) = {23};\n')
    fid.write('%s(25) = {-7,4,9};\n' % curvename)
    fid.write('Surface(26) = {25};\n')
    fid.write('%s(27) = {-4,12,-6};\n' % curvename)
    fid.write('Surface(28) = {27};\n')
    fid.write('Surface Loop(29) = {28,26,16,14,20,24,22,18};\n')
    fid.write('Volume(30) = {29};\n')
    fid.write('Physical Surface(1) = {28,26,16,14,20,24,22,18};\n')
    fid.write('Physical Volume(2) = 30;\n')
    fid.close()
    # }}}

    if options.exist('refine'):
        meshini = options.getfieldvalue('refine')
        metric = options.getfieldvalue('refinemetric')

        # Create .pos file with existing mesh and refining metric:  {{{
        fid = open('sphere.pos', 'w')

        fid.write('View "background mesh" {\n')
        for i in range(meshini.numberofelements):
            fid.write('ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n'
                      % (meshini.x[meshini.elements[i, 0] - 1], meshini.y[meshini.elements[i, 0] - 1], meshini.z[meshini.elements[i, 0] - 1],
                         meshini.x[meshini.elements[i, 1] - 1], meshini.y[meshini.elements[i, 1] - 1], meshini.z[meshini.elements[i, 1] - 1],
                         meshini.x[meshini.elements[i, 2] - 1], meshini.y[meshini.elements[i, 2] - 1], meshini.z[meshini.elements[i, 2] - 1],
                         metric[meshini.elements[i, 0] - 1], metric[meshini.elements[i, 1] - 1], metric[meshini.elements[i, 2] - 1]))
        fid.write('};\n')
        fid.close()
        # }}}

    # Call gmsh
    #
    # NOTE: The default format in Gmsh 3 is "msh2". Rather than conditionally
    #       modifying our parsing scheme for Gmsh 4, for now, we simply set the
    #       "-format" option.
    #
    if options.exist('refine'):
        subprocess.call('gmsh -2 sphere.geo -bgm sphere.pos', shell=True)
    else:
        subprocess.call('gmsh -2 sphere.geo', shell=True)

    # Import mesh  {{{
    fid = open('sphere.msh', 'r')

    # Get mesh format
    A = fid.readline().strip()
    if A != '$MeshFormat':
        raise RuntimeError(['Expecting $MeshFormat (', A, ')'])

    A = fid.readline().split()
    A = fid.readline().strip()
    if A != '$EndMeshFormat':
        raise RuntimeError(['Expecting $EndMeshFormat (', A, ')'])

    # Nodes
    A = fid.readline().strip()
    if A != '$Nodes':
        raise RuntimeError(['Expecting $Nodes (', A, ')'])

    mesh.numberofvertices = int(fid.readline().strip())
    mesh.x = np.empty(mesh.numberofvertices)
    mesh.y = np.empty(mesh.numberofvertices)
    mesh.z = np.empty(mesh.numberofvertices)
    for i in range(mesh.numberofvertices):
        A = fid.readline().split()
        mesh.x[i] = float(A[1])
        mesh.y[i] = float(A[2])
        mesh.z[i] = float(A[3])

    A = fid.readline().strip()
    if A != '$EndNodes':
        raise RuntimeError(['Expecting $EndNodes (', A, ')'])

    # Elements
    A = fid.readline().strip()
    if A != '$Elements':
        raise RuntimeError(['Expecting $Elements (', A, ')'])
    mesh.numberofelements = int(fid.readline().strip())
    mesh.elements = np.zeros([mesh.numberofelements, 3])
    for i in range(mesh.numberofelements):
        A = fid.readline().split()
        mesh.elements[i] = [int(A[5]), int(A[6]), int(A[7])]
    mesh.elements = mesh.elements.astype(int)
    A = fid.readline().strip()
    if A != '$EndElements':
        raise RuntimeError(['Expecting $EndElements (', A, ')'])
    fid.close()
    # }}}

    # A little technicality here. The mesh is not exactly on the sphere. We
    # create lat,long coordinates, and reproject onto an exact sphere.
    mesh.r = np.sqrt(mesh.x ** 2 + mesh.y ** 2 + mesh.z ** 2)

    # Make sure we don't have south and north pole
    pos = np.where(np.logical_and.reduce((mesh.x == 0, mesh.y == 0)))[0]
    mesh.lat = asind(mesh.z / mesh.r)
    mesh.long = atan2d(mesh.y, mesh.x)
    pos = np.where(mesh.lat == 90)[0]
    mesh.lat[pos] = 90 - 0.01
    pos = np.where(mesh.lat == -90)[0]
    mesh.lat[pos] = -90 + 0.01

    mesh.r = radius * np.ones((mesh.numberofvertices, ))
    mesh.x = radius * cosd(mesh.lat) * cosd(mesh.long)
    mesh.y = radius * cosd(mesh.lat) * sind(mesh.long)
    mesh.z = radius * sind(mesh.lat)

    # Erase files
    subprocess.call('rm -rf sphere.geo sphere.msh sphere.pos', shell=True)

    # Return mesh
    return mesh
