#Test Name: RoundSheetShelfGLMigrationSSA2d
import numpy as np
from model import *
from socket import gethostname
from roundmesh import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


radius = 1.e6
shelfextent = 2.e5

md = roundmesh(model(), radius, 50000.)
#fix center node to 0, 0
rad = np.sqrt(md.mesh.x**2 + md.mesh.y**2)
pos = np.argmin(rad)
md.mesh.x[pos] = 0.
md.mesh.y[pos] = 0.  #the closest node to the center is changed to be exactly at the center
xelem = np.mean(md.mesh.x[md.mesh.elements - 1], axis=1)
yelem = np.mean(md.mesh.y[md.mesh.elements - 1], axis=1)
rad = np.sqrt(xelem**2 + yelem**2)
flags = np.zeros(md.mesh.numberofelements)
pos = np.nonzero(rad >= (radius - shelfextent))
flags[pos] = 1
md = setmask(md, flags, '')
md = parameterize(md, '../Par/RoundSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

md.transient.isthermal = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isstressbalance = True
md.transient.isgroundingline = True

#test different grounding line dynamics.
md.groundingline.migration = 'AggressiveMigration'
md = solve(md, 'Transient')
element_on_iceshelf_agressive = md.results.TransientSolution[0].MaskOceanLevelset
vel_agressive = md.results.TransientSolution[0].Vel

md.groundingline.migration = 'SoftMigration'
md = solve(md, 'Transient')
element_on_iceshelf_soft = md.results.TransientSolution[0].MaskOceanLevelset
vel_soft = md.results.TransientSolution[0].Vel

md.mask.ocean_levelset = md.geometry.thickness + md.materials.rho_water / md.materials.rho_ice * md.geometry.bed
md.groundingline.migration = 'SubelementMigration'
md.groundingline.friction_interpolation = 'SubelementFriction1'
md = solve(md, 'Transient')
element_on_iceshelf_subelement = md.results.TransientSolution[0].MaskOceanLevelset
vel_subelement = md.results.TransientSolution[0].Vel

md.groundingline.friction_interpolation = 'SubelementFriction2'
md = solve(md, 'Transient')
element_on_iceshelf_subelement2 = md.results.TransientSolution[0].MaskOceanLevelset
vel_subelement2 = md.results.TransientSolution[0].Vel

#Fields and tolerances to track changes
field_names = ['ElementOnIceShelfAggressive', 'VelAggressive', 'ElementOnIceShelfSoft', 'VelSoft', 'ElementOnIceShelfSubelement', 'VelSubelement', 'ElementOnIceShelfSubelement2', 'VelSubelement2']
field_tolerances = [1e-13, 2e-13, 1e-13, 2e-13, 1e-13, 2e-13, 1e-13, 2e-13]
field_values = [element_on_iceshelf_agressive, vel_agressive, element_on_iceshelf_soft, vel_soft, element_on_iceshelf_subelement, vel_subelement, element_on_iceshelf_subelement2, vel_subelement2]
