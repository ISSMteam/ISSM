#Test Name: Esa2Dsurface
#Elastostatic adjustment for an elemental ice unloading
from socket import gethostname

import numpy as np

from lovenumbers import *
from model import *
from paterson import *
from roundmesh import *
from solve import *


#mesh earth:
md = model()
md = roundmesh(md, 50000, 2000)  # radius and element size (meters)

#define load
md.esa.deltathickness = np.zeros((md.mesh.numberofelements, ))
disc_radius = 20  # km
index = md.mesh.elements
x_element = np.mean(md.mesh.x[index - 1], 1)
y_element = np.mean(md.mesh.y[index - 1], 1)
rad_dist = np.sqrt(x_element**2 + y_element**2) / 1000  # radial distance in km
md.esa.deltathickness[np.where(rad_dist <= disc_radius)] = -1  # 1 m water withdrawal

#love numbers:
md.solidearth.lovenumbers = lovenumbers('maxdeg', 10000)

#mask:  {{{
#make sure wherever there is an ice load, that the mask is set to ice:
md.mask.ice_levelset = np.ones((md.mesh.numberofvertices, ))
pos = np.where(md.esa.deltathickness)
md.mask.ice_levelset[md.mesh.elements[pos, :] - 1] = -1

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
md.miscellaneous.name = 'test2110'
# }}}

md.esa.degacc = 0.01

#solve esa
md.esa.requested_outputs = ['EsaUmotion', 'EsaXmotion', 'EsaYmotion',
                            'EsaStrainratexx', 'EsaStrainratexy', 'EsaStrainrateyy', 'EsaRotationrate']
md.cluster = generic('name', gethostname(), 'np', 3)
md.verbose = verbose('111111111')
md = solve(md, 'Esa')

#Fields and tolerances to track changes
field_names = ['EsaUmotion', 'EsaXmotion', 'EsaYmotion',
               'EsaStrainratexx', 'EsaStrainratexy', 'EsaStrainrateyy', 'EsaRotationrate']
field_tolerances = [1e-13, 2e-12, 2e-12, 9e-12, 8e-12, 8e-12, 3e-11]
field_values = [md.results.EsaSolution.EsaUmotion,
                md.results.EsaSolution.EsaXmotion,
                md.results.EsaSolution.EsaYmotion,
                md.results.EsaSolution.EsaStrainratexx,
                md.results.EsaSolution.EsaStrainratexy,
                md.results.EsaSolution.EsaStrainrateyy,
                md.results.EsaSolution.EsaRotationrate]
