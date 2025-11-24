#Test Name: ValleyGlacierLevelsetEnthalpyHO3d
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 50000)
md = setmask(md, '', '')
md = parameterize(md, '../Par/ValleyGlacierShelf.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#Thermal model
pos_surf = np.where(md.mesh.vertexonsurface)[0]
md.thermal.spctemperature[pos_surf] = md.initialization.temperature[pos_surf]
md.thermal.isenthalpy = True
md.thermal.isdynamicbasalspc = True

#Transient
md.transient.isstressbalance = True
md.transient.ismovingfront = True
md.transient.ismasstransport = True
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = True
md.groundingline.melt_interpolation = 'SubelementMelt1'

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Thickness1', 'Surface1', 'MaskIceLevelset1', 'Enthalpy1', 'Watercolumn1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Thickness2', 'Surface2', 'MaskIceLevelset2', 'Enthalpy2', 'Watercolumn2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Thickness3', 'Surface3', 'MaskIceLevelset3', 'Enthalpy3', 'Watercolumn3']
field_tolerances = [2e-12, 2e-12, 2e-12, 3e-13, 3e-13, 2e-13, 2e-13, 1e-13, 1e-13,
                    4e-12, 6e-12, 4e-12, 3e-13, 8e-13, 2e-13, 2e-13, 1e-13, 1e-13,
                    6e-12, 2e-11, 6e-12, 8e-13, 9e-13, 3e-13, 2e-13, 1e-13, 1e-13]
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
