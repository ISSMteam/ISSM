#Test Name: EarthSlc
import numpy as np
import pickle
#from gmshplanet import *
from gmtmask import *
from lovenumbers import *
from materials import *
from model import *
from paterson import *
from solve import *


# Mesh earth
#
# NOTE: In MATLAB, we currently use cached mesh to account for differences in
# mesh generated under Linux versus under macOS
#
md = model()

# Generate and save mesh (need to uncomment import of gmshplanet as well)
# md.mesh = gmshplanet('radius', 6.371012 * 1e3, 'resolution', 700.) # 700 km resolution mesh
# with open('../Data/SlcTestMesh.pkl', 'wb') as slc_test_mesh_file:
#     pickle.dump(md.mesh, slc_test_mesh_file)

# Load precomputed mesh
with open('../Data/SlcTestMesh.pkl', 'rb') as slc_test_mesh_file:
    md.mesh = pickle.load(slc_test_mesh_file, encoding='latin1')

# Geometry for the bed, arbitrary thickness of 100
md.geometry.bed = np.zeros((md.mesh.numberofvertices, ))
md.geometry.base = md.geometry.bed
md.geometry.thickness = 100 * np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = md.geometry.bed + md.geometry.thickness

# Solidearth loading #{{{
md.masstransport.spcthickness = np.append(md.geometry.thickness, 0)
md.smb.mass_balance = np.zeros((md.mesh.numberofvertices, 1))

# Antarctica
xe = md.mesh.x[md.mesh.elements - 1].sum(axis=1) / 3
ye = md.mesh.y[md.mesh.elements - 1].sum(axis=1) / 3
ze = md.mesh.z[md.mesh.elements - 1].sum(axis=1) / 3
re = pow((pow(xe, 2) + pow(ye, 2) + pow(ze, 2)), 0.5)

late = asind(ze / re)
longe = atan2d(ye, xe)
pos = np.where(late < -80)[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 100
posant = pos

# Greenland
pos = np.where(np.logical_and.reduce((late > 60, late < 90, longe > -75, longe < -15)))[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 100
posgre = pos

# Elastic loading from love numbers:
md.solidearth.lovenumbers = lovenumbers('maxdeg', 1000)
#}}}

# Mask: {{{
mask = gmtmask(md.mesh.lat, md.mesh.long)
oceanmask = -1 * np.ones((md.mesh.numberofvertices, ))
pos = np.where(mask == 0)[0]
oceanmask[pos] = 1

icemask = np.ones((md.mesh.numberofvertices, ))
# NOTE: Need to be careful here: when using multidimensional array elements to
# address a one-dimensional array in MATLAB, only first column is used
#
icemask[md.mesh.elements[posant][:, 0] - 1] = -1
icemask[md.mesh.elements[posgre][:, 0] - 1] = -1

md.mask.ice_levelset = icemask
md.mask.ocean_levelset = oceanmask
#}}}

# Time stepping {{{
md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 1
#}}}

# Masstransport
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vx = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vy = np.zeros((md.mesh.numberofvertices, ))
md.initialization.sealevel = np.zeros((md.mesh.numberofvertices, ))
md.initialization.str = 0

# Materials
md.materials = materials('hydro')

# Miscellaneous
md.miscellaneous.name = 'test2002'

# Solution parameters
md.cluster.np = 3
md.solidearth.settings.reltol = np.nan
md.solidearth.settings.abstol = 1e-3
md.solidearth.settings.sealevelloading = 1
md.solidearth.settings.isgrd = 1
md.solidearth.settings.ocean_area_scaling = 0
md.solidearth.settings.grdmodel = 1

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1
md.solidearth.requested_outputs = ['Sealevel', 'Bed', 'SealevelBarystaticIceLoad', 'SealevelBarystaticIceArea', 'SealevelBarystaticIceWeights']

# Max number of iterations reverted back to 10 (i.e., the original default value)
md.solidearth.settings.maxiter = 10

# Eustatic run
md.solidearth.settings.selfattraction = 0
md.solidearth.settings.elastic = 0
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0

md = solve(md, 'Transient')
Seustatic = md.results.TransientSolution.Sealevel
Beustatic = md.results.TransientSolution.Bed

# Eustatic + selfattraction run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 0
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md = solve(md, 'tr')
Sselfattraction = md.results.TransientSolution.Sealevel
Bselfattraction = md.results.TransientSolution.Bed

# Eustatic + selfattraction + elastic run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md = solve(md, 'tr')
Selastic = md.results.TransientSolution.Sealevel
Belastic = md.results.TransientSolution.Bed

# Eustatic + selfattraction + elastic + rotation run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 1
md.solidearth.settings.viscous = 0
md = solve(md, 'tr')
Srotation = md.results.TransientSolution.Sealevel
Brotation = md.results.TransientSolution.Bed

# Fields and tolerances to track changes
field_names = ['Seustatic', 'Sselfattraction', 'Selastic', 'Srotation', 'Beustatic', 'Bselfattraction', 'Belastic', 'Brotation']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [Seustatic, Sselfattraction, Selastic, Srotation, Beustatic, Bselfattraction, Belastic, Brotation]
