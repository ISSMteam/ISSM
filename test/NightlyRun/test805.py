#Test Name: ValleyGlacierLevelsetEnthCalvingHO3d
import numpy as np
from socket import gethostname
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *


md = triangle(model(), '../Exp/Square.exp', 50000)
md = setmask(md, '', '')
md = parameterize(md, '../Par/ValleyGlacierShelf.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#Thermal model
pos_surf = np.nonzero(md.mesh.vertexonsurface)[0]
md.thermal.spctemperature[pos_surf] = md.initialization.temperature[pos_surf]
md.thermal.isenthalpy = True
md.thermal.isdynamicbasalspc = True

#Transient
md.transient.isstressbalance = True
md.transient.ismasstransport = True
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = True
md.transient.ismovingfront = True

md.calving.calvingrate = 1000. * np.ones((md.mesh.numberofvertices))
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices))
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.levelset.stabilization = 2
md.levelset.migration_max = 1e10

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1', 'Enthalpy1', 'Watercolumn1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2', 'Enthalpy2', 'Watercolumn2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3', 'Enthalpy3', 'Watercolumn3']
field_tolerances = [1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11,
                    1e-9, 1e-9, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-9, 1e-10,
                    1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].MaskIceLevelset,
                md.results.TransientSolution[0].Enthalpy,
                md.results.TransientSolution[0].Watercolumn,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].MaskIceLevelset,
                md.results.TransientSolution[1].Enthalpy,
                md.results.TransientSolution[1].Watercolumn,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].MaskIceLevelset,
                md.results.TransientSolution[2].Enthalpy,
                md.results.TransientSolution[2].Watercolumn]
