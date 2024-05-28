#Test Name: FlowbandFSsheetshelf
import numpy as np
#from scipy.interpolate import interp1d
from model import *
from setflowequation import *
from solve import *
from NowickiProfile import *
from bamgflowband import *
from paterson import *

#mesh parameters
x = np.arange(-5, 5.5, .5).T
[b, h, sea] = NowickiProfile(x)
x = x * 10**3
h = h * 10**3
b = (b - sea) * 10**3

#mesh domain
md = bamgflowband(model(), x, b + h, b, 'hmax', 150.)

#parameterize
md.geometry.surface = np.interp(md.mesh.x, x, b + h)
md.geometry.base = np.interp(md.mesh.x, x, b)
md.geometry.thickness = md.geometry.surface - md.geometry.base

#mask
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices, ))
md.mask.ice_levelset[np.where(md.mesh.vertexflags(2))] = 0
md.mask.ocean_levelset = -0.5 * np.ones((md.mesh.numberofvertices))
md.mask.ocean_levelset[np.where(md.mesh.x < 0)] = 0.5

#materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3 * np.ones((md.mesh.numberofelements, ))

#damage
md.damage.D = np.zeros((md.mesh.numberofvertices, ))
md.damage.spcdamage = float('NaN') * np.ones((md.mesh.numberofvertices, ))

#friciton
md.friction.coefficient = np.zeros((md.mesh.numberofvertices, ))
md.friction.coefficient[np.where(md.mesh.vertexflags(1))] = 20
md.friction.p = np.ones((md.mesh.numberofelements, ))
md.friction.q = np.ones((md.mesh.numberofelements, ))

#boundary conditions
md.stressbalance.spcvx = float('NaN') * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvy = float('NaN') * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvz = float('NaN') * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.referential = float('NaN') * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = np.zeros((md.mesh.numberofvertices, 3))
md.stressbalance.spcvx[np.where(md.mesh.vertexflags(4))] = 800.
md.stressbalance.spcvy[np.where(md.mesh.vertexflags(4))] = 0.
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))

#misc
md = setflowequation(md, 'FS', 'all')
md.stressbalance.abstol = float('NaN')
md.stressbalance.FSreconditioning = 1
md.stressbalance.maxiter = 20
md.flowequation.augmented_lagrangian_r = 10000.
md.flowequation.augmented_lagrangian_rhop = 10000.
md.initialization.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.y)
md.miscellaneous.name = 'test702'
md.groundingline.migration = 'None'
md.cluster = generic('np', 2)

#Fields and tolerances to track changes
field_names = []
field_tolerances = []
field_values = []
for i in ['MINI', 'MINIcondensed', 'TaylorHood', 'XTaylorHood', 'LATaylorHood']:
    print(' ')
    print(' == == == Testing ' + i + ' Full - Stokes Finite element == == = ')
    md.flowequation.fe_FS = i
    md = solve(md, 'Stressbalance')
    field_names.extend(['Vx' + i, 'Vy' + i, 'Vel' + i, 'Pressure' + i])
    field_tolerances.extend([8e-5, 8e-5, 8e-5, 1e-08])
    field_values.extend([md.results.StressbalanceSolution.Vx,
                         md.results.StressbalanceSolution.Vy,
                         md.results.StressbalanceSolution.Vel,
                         md.results.StressbalanceSolution.Pressure])
