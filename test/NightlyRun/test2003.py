#Test Name: EarthSlc_rotationalFeedback
import numpy as np
import pickle
from socket import gethostname
from gmtmask import *
from lovenumbers import *
from MatlabFuncs import *
from model import *
from paterson import *
from solve import *


md = model()

# Load precomputed mesh
with open('../Data/SlcTestMesh.pkl', 'rb') as slc_test_mesh_file:
    md.mesh = pickle.load(slc_test_mesh_file, encoding='latin1')

# Geometry for the bed, arbitrary thickness of 100
md.geometry.bed = -1 * np.ones((md.mesh.numberofvertices, ))
md.geometry.base = md.geometry.bed
md.geometry.thickness = 1000 * np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = md.geometry.bed + md.geometry.thickness

# Parameterize slc solution:
#solidearth loading:  {{{
md.masstransport.spcthickness = np.append(md.geometry.thickness, 0)
md.smb.mass_balance = np.zeros(md.mesh.numberofvertices)

xe = md.mesh.x[md.mesh.elements - 1].sum(axis=1) / 3
ye = md.mesh.y[md.mesh.elements - 1].sum(axis=1) / 3
ze = md.mesh.z[md.mesh.elements - 1].sum(axis=1) / 3
re = pow((pow(xe, 2) + pow(ye, 2) + pow(ze, 2)), 0.5)

late = asind(ze / re)
longe = atan2d(ye, xe)

# Greenland
pos = np.where(np.logical_and.reduce((late > 60, late < 90, longe > -75, longe < -15)))[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 1000
posice = pos

# Elastic loading from love numbers
md.solidearth.lovenumbers = lovenumbers('maxdeg', 100)
#}}}

# Mask: {{{
mask = gmtmask(md.mesh.lat, md.mesh.long)
icemask = np.ones(md.mesh.numberofvertices)
icemask[md.mesh.elements[posice, :] - 1] = -1
md.mask.ice_levelset = icemask
oceanmask = -1 * np.ones(md.mesh.numberofvertices)
pos = np.where(mask == 0)[0]
oceanmask[pos] = 1
md.mask.ocean_levelset = oceanmask

# Use model representation of ocean area (not the true area)
md.solidearth.settings.ocean_area_scaling = 0

# Materials
md.initialization.temperature = 273.25 * np.ones(md.mesh.numberofvertices)
md.initialization.sealevel = np.zeros(md.mesh.numberofvertices)
md.initialization.str = 0

md.basalforcings.groundedice_melting_rate = np.zeros(md.mesh.numberofvertices)
md.basalforcings.floatingice_melting_rate = np.zeros(md.mesh.numberofvertices)
md.initialization.vx = np.zeros(md.mesh.numberofvertices)
md.initialization.vy = np.zeros(md.mesh.numberofvertices)

# Miscellaneous
md.miscellaneous.name = 'test2003'

# Solution parameters
md.solidearth.settings.reltol = np.nan
md.solidearth.settings.abstol = 1e-3
md.solidearth.settings.sealevelloading = 0
md.solidearth.settings.grdocean = 0
md.solidearth.settings.isgrd = 1
md.solidearth.settings.ocean_area_scaling = 0
md.solidearth.settings.grdmodel = 1
md.solidearth.settings.horiz = 1
md.solidearth.requested_outputs = ['Sealevel', 'Bed', 'BedEast', 'BedNorth']

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1

md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 1

# Eustatic + selfattraction + elastic run:
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md.cluster = generic('name', gethostname(), 'np', 3)
#md.verbose = verbose('111111111')
md = solve(md, 'Transient')
SnoRotation = md.results.TransientSolution.Sealevel
BUnoRotation = md.results.TransientSolution.Bed
BEnoRotation = md.results.TransientSolution.BedEast
BNnoRotation = md.results.TransientSolution.BedNorth

# Eustatic + selfattraction + elastic + rotation run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 1
md.solidearth.settings.viscous = 0
md.cluster = generic('name', gethostname(), 'np', 3)
#md.verbose = verbose('111111111')
md = solve(md, 'Transient')
SRotation = md.results.TransientSolution.Sealevel
BURotation = md.results.TransientSolution.Bed
BERotation = md.results.TransientSolution.BedEast
BNRotation = md.results.TransientSolution.BedNorth

# Fields and tolerances to track changes
field_names = ['Sealevel', 'Uplift', 'NorthDisplacement', 'EastDisplacement']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [SRotation - SnoRotation, BURotation - BUnoRotation, BNRotation - BNnoRotation,BERotation - BEnoRotation]
