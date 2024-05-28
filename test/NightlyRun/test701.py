#Test Name: FlowbandFSshelf
import numpy as np
from model import *
from solve import *
from setflowequation import *
from bamgflowband import *
from paterson import *

x = np.arange(1, 3001, 100).T
h = np.linspace(1000, 300, np.size(x)).T
b = -917. / 1023. * h

md = bamgflowband(model(), x, b + h, b, 'hmax', 80.)

#Geometry  #interp1d returns a function to be called on md.mesh.x
md.geometry.surface = np.interp(md.mesh.x, x, b + h)
md.geometry.base = np.interp(md.mesh.x, x, b)
md.geometry.thickness = md.geometry.surface - md.geometry.base

#mask
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices))
md.mask.ice_levelset[np.where(md.mesh.vertexflags(2))] = 0.
md.mask.ocean_levelset = np.zeros((md.mesh.numberofvertices)) - 0.5

#materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#friction
md.friction.coefficient = np.zeros((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mesh.vertexflags(1))] = 20.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

#Boundary conditions
md.stressbalance.referential = np.nan * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = 0. * np.ones((md.mesh.numberofvertices, 3))
md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvx[np.where(md.mesh.vertexflags(4))] = 0.
md.stressbalance.spcvy[np.where(md.mesh.vertexflags(4))] = 0.
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))

#Misc
md = setflowequation(md, 'FS', 'all')
md.stressbalance.abstol = np.nan
#md.stressbalance.reltol = 10**-16
md.stressbalance.FSreconditioning = 1.
md.stressbalance.maxiter = 20
md.flowequation.augmented_lagrangian_r = 10000.
md.miscellaneous.name = 'test701'
md.verbose = verbose('convergence', True)
md.cluster = generic('np', 2)
md.groundingline.migration = 'None'

#Go solve
field_names = []
field_tolerances = []
field_values = []
#md.initialization.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.y)
for i in ['MINI', 'MINIcondensed', 'TaylorHood', 'LATaylorHood', 'CrouzeixRaviart', 'LACrouzeixRaviart']:
    print(' ')
    print('=====Testing ' + i + ' Full-Stokes Finite element=====')
    md.flowequation.fe_FS = i
    md = solve(md, 'Stressbalance')
    field_names = field_names + ['Vx' + i, 'Vy' + i, 'Vel' + i, 'Pressure' + i]
    field_tolerances = field_tolerances + [9e-5, 9e-5, 9e-5, 1e-10]
    field_values = field_values + [md.results.StressbalanceSolution.Vx,
                                   md.results.StressbalanceSolution.Vy,
                                   md.results.StressbalanceSolution.Vel,
                                   md.results.StressbalanceSolution.Pressure]
