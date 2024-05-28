#Test Name: MomentOfIntertia
import numpy as np
import pickle
from socket import gethostname
from gmtmask import *
from lovenumbers import *
from model import *
from paterson import *
from solve import *


# Mesh earth
md = model()

# Load precomputed mesh
with open('../Data/SlcTestMesh.pkl', 'rb') as slc_test_mesh_file:
    md.mesh = pickle.load(slc_test_mesh_file, encoding='latin1')

# Geometry for the bed, arbitrary thickness of 100
md.geometry.bed = -1 * np.ones((md.mesh.numberofvertices, ))
md.geometry.base = md.geometry.bed
md.geometry.thickness = 100 * np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = md.geometry.bed + md.geometry.thickness

# Parameterize SLC solution
# Solidearth loading  {{{
md.masstransport.spcthickness = np.append(md.geometry.thickness, 0)
md.smb.mass_balance = np.zeros((md.mesh.numberofvertices, ))

xe = md.mesh.x[md.mesh.elements - 1].sum(axis=1) / 3
ye = md.mesh.y[md.mesh.elements - 1].sum(axis=1) / 3
ze = md.mesh.z[md.mesh.elements - 1].sum(axis=1) / 3
re = pow((pow(xe, 2) + pow(ye, 2) + pow(ze, 2)), 0.5)

late = asind(ze / re)
longe = atan2d(ye, xe)
# Greenland
pos = np.where(np.logical_and.reduce((late > 60, late < 90, longe > -75, longe < -15)))[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 100
posice = pos

# Elastic loading from love numbers
md.solidearth.lovenumbers = lovenumbers('maxdeg', 100)

#}}}
#mask:  {{{
mask = gmtmask(md.mesh.lat, md.mesh.long)
icemask = np.ones((md.mesh.numberofvertices, 1))
icemask[md.mesh.elements[posice, :] - 1] = -0.5

oceanmask = -1 * np.ones((md.mesh.numberofvertices, 1))
pos = np.where(mask == 0)[0]
oceanmask[pos] = 1
icemask[np.logical_not(pos).astype(int)] = 1

md.mask.ice_levelset = icemask
md.mask.ocean_levelset = oceanmask

# Use model representation of ocean area (not the true area)
md.solidearth.settings.ocean_area_scaling = 0

# Materials
md.initialization.temperature = 273.25 * np.ones((md.mesh.numberofvertices, 1))
md.initialization.sealevel = np.zeros((md.mesh.numberofvertices, 1))
md.initialization.str = 0

md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vx = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vy = np.zeros((md.mesh.numberofvertices, ))

# Miscellaneous
md.miscellaneous.name = 'test2010'

# Solution parameters
md.solidearth.settings.reltol = np.nan
md.solidearth.settings.abstol = 1e-3
md.solidearth.settings.sealevelloading = 0
md.solidearth.settings.grdocean = 1
md.solidearth.settings.isgrd = 1
md.solidearth.settings.ocean_area_scaling = 0
md.solidearth.settings.grdmodel = 1
md.solidearth.settings.horiz = 1
md.solidearth.requested_outputs = [
    'Sealevel',
    'SealevelBarystaticIceArea',
    'SealevelBarystaticIceLoad',
    'SealevelBarystaticIceMask',
    'SealevelBarystaticIceLatbar',
    'SealevelBarystaticIceLongbar'
]

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1

md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 1

# Eustatic + selfattraction + elastic + rotation run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 1
md.solidearth.settings.viscous = 0
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

moi_p = md.solidearth.rotational.polarmoi
moi_e = md.solidearth.rotational.equatorialmoi
tide_love_k2 = md.solidearth.lovenumbers.tk[2]
load_love_k2 = md.solidearth.lovenumbers.k[2]
tide_love_k2secular = md.solidearth.lovenumbers.tk2secular
# uncomment following 2 lines for
eus = md.results.TransientSolution.Bslc
slc = md.results.TransientSolution.Sealevel
moixz = md.results.TransientSolution.SealevelchangePolarMotionX / (1 / (1 - tide_love_k2 / tide_love_k2secular) * (1 + load_love_k2) / (moi_p - moi_e))
moiyz = md.results.TransientSolution.SealevelchangePolarMotionY / (1 / (1 - tide_love_k2 / tide_love_k2secular) * (1 + load_love_k2) / (moi_p - moi_e))
moizz = md.results.TransientSolution.SealevelchangePolarMotionZ / ( -(1 + load_love_k2) / moi_p)

areaice = md.results.TransientSolution.SealevelBarystaticIceArea
areaice[np.isnan(areaice)] = 0
loadice = md.results.TransientSolution.SealevelBarystaticIceLoad
rad_e = md.solidearth.planetradius

lat = md.results.TransientSolution.SealevelBarystaticIceLatbar * np.pi / 180
lon = md.results.TransientSolution.SealevelBarystaticIceLongbar * np.pi / 180
moi_xz = sum(-loadice * areaice * pow(rad_e, 2) * np.sin(lat) * np.cos(lat) * np.cos(lon))
moi_yz = sum(-loadice * areaice * pow(rad_e, 2) * np.sin(lat) * np.cos(lat) * np.sin(lon))
moi_zz = sum(-loadice * areaice * pow(rad_e, 2) * (-1.0 / 3.0 + np.sin(lat) ** 2))
theoretical_value_check = [moixz / moi_xz, moiyz / moi_yz, moizz / moi_zz] # Should yield [1.0, 1.0, 1.0]
# }}}

#Fields and tolerances to track changes
field_names = ['eus', 'slc', 'moixz', 'moiyz', 'moizz']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [eus, slc, moixz, moiyz, moizz]
