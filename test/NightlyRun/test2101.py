#Test Name: EarthEsa
#Elastostatic adjustment for an elemental ice unloading
import numpy as np
import pickle
from socket import gethostname
from gmtmask import *
from lovenumbers import *
from model import *
from paterson import *
from solve import *


#mesh earth:
md = model()

# Load precomputed mesh
with open('../Data/SlcTestMesh.pkl', 'rb') as slc_test_mesh_file:
    md.mesh = pickle.load(slc_test_mesh_file, encoding='latin1')

#define load
md.esa.deltathickness = np.zeros((md.mesh.numberofelements, ))
pos = 449
md.esa.deltathickness[pos] = -100  # this is the only "icy" element

#love numbers:
md.solidearth.lovenumbers = lovenumbers('maxdeg', 10000)

#mask:  {{{
md.mask.ocean_levelset = gmtmask(md.mesh.lat, md.mesh.long)

#make sure wherever there is an ice load, that the mask is set to ice:
md.mask.ice_levelset = np.ones((md.mesh.numberofvertices, ))
pos = np.where(md.esa.deltathickness)
md.mask.ice_levelset[md.mesh.elements[pos, :]] = -1

#is ice grounded?
md.mask.ocean_levelset = -np.ones((md.mesh.numberofvertices, ))
pos = np.where(md.mask.ice_levelset <= 0)
md.mask.ocean_levelset[pos] = 1
# }}}
#geometry:  {{{
di = md.materials.rho_ice / md.materials.rho_water
md.geometry.thickness = np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = (1 - di) * np.zeros((md.mesh.numberofvertices, ))
md.geometry.base = md.geometry.surface - md.geometry.thickness
md.geometry.bed = md.geometry.base
# }}}
#materials:  {{{
md.initialization.temperature = 273.25 * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements, ))
# }}}
#Miscellaneous: {{{
md.miscellaneous.name = 'test2101'
# }}}

#solve esa
md.esa.requested_outputs = ['EsaUmotion', 'EsaNmotion', 'EsaEmotion']
md.cluster = generic('name', gethostname(), 'np', 3)
md.verbose = verbose('111111111')
md = solve(md, 'Esa')

#Fields and tolerances to track changes
field_names = ['EsaUmotion', 'EsaNmotion', 'EsaEmotion']
field_tolerances = [1e-13, 1e-13, 2e-13]
field_values = [md.results.EsaSolution.EsaUmotion,
                md.results.EsaSolution.EsaNmotion,
                md.results.EsaSolution.EsaEmotion]
