#Test Name: ThermalEnthBasalcondsTrans
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 300000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareThermal.py')

h = 100.
md.geometry.thickness = h * np.ones((md.mesh.numberofvertices, ))
md.geometry.base = -h * np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = md.geometry.base + md.geometry.thickness

md.extrude(41, 2.)
md = setflowequation(md, 'HO', 'all')
md.thermal.isenthalpy = True
md.thermal.isdynamicbasalspc = True

#Basal forcing
Ts = 273.15 - 3.
Tb = 273.15 - 1.
Tsw = Tb
qgeo = md.materials.thermalconductivity / max(md.geometry.thickness) * (Tb - Ts)  #qgeo = kappa * (Tb - Ts) / H
md.basalforcings.geothermalflux[np.where(md.mesh.vertexonbase)] = qgeo
md.initialization.temperature = qgeo / md.materials.thermalconductivity * (md.geometry.surface - md.mesh.z) + Ts

#Surface forcing
pos = np.where(md.mesh.vertexonsurface)
SPC_cold = float('NaN') * np.ones((md.mesh.numberofvertices, ))
SPC_warm = float('NaN') * np.ones((md.mesh.numberofvertices, ))
SPC_cold[pos] = Ts
SPC_warm[pos] = Tsw
md.thermal.spctemperature = SPC_cold
md.timestepping.time_step = 5.
t0 = 0.
t1 = 10.
t2 = 100.
md.timestepping.final_time = 400.
md.thermal.spctemperature = np.array([SPC_cold, SPC_cold, SPC_warm, SPC_warm, SPC_cold]).T
md.thermal.spctemperature = np.vstack((md.thermal.spctemperature, [t0, t1 - 1, t1, t2 - 1, t2]))
#print np.shape(md.thermal.spctemperature)
#print md.thermal.spctemperature

#Additional settings
md.transient.ismasstransport = False
md.transient.isstressbalance = False
md.transient.issmb = True
md.transient.isthermal = True
md.thermal.stabilization = False

#Go solve
md.verbose = verbose(0)
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Enthalpy1', 'Temperature1', 'Waterfraction1', 'BasalMeltingRate1', 'Watercolumn1',
               'Enthalpy2', 'Temperature2', 'Waterfraction2', 'BasalMeltingRate2', 'Watercolumn2',
               'Enthalpy3', 'Temperature3', 'Waterfraction3', 'BasalMeltingRate3', 'Watercolumn3',
               'Enthalpy4', 'Temperature4', 'Waterfraction4', 'BasalMeltingRate4', 'Watercolumn4']
field_tolerances = [1.e-10, 1.e-10, 1.e-10, 1.e-9, 1.e-10,
                    1.e-10, 1.e-10, 1.e-10, 2.e-9, 2.e-10,
                    1.e-10, 1.e-10, 1.e-10, 2.e-9, 1.e-10,
                    1.e-10, 1.e-10, 1.e-10, 2.e-9, 1.e-10]
i1 = 0
i2 = int(np.ceil(t2 / md.timestepping.time_step) + 2) - 1
i3 = int(np.ceil(md.timestepping.final_time / (2. * md.timestepping.time_step))) - 1
i4 = len(md.results.TransientSolution) - 1
field_values = [md.results.TransientSolution[i1].Enthalpy,
                md.results.TransientSolution[i1].Temperature,
                md.results.TransientSolution[i1].Waterfraction,
                md.results.TransientSolution[i1].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[i1].Watercolumn,
                md.results.TransientSolution[i2].Enthalpy,
                md.results.TransientSolution[i2].Temperature,
                md.results.TransientSolution[i2].Waterfraction,
                md.results.TransientSolution[i2].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[i2].Watercolumn,
                md.results.TransientSolution[i3].Enthalpy,
                md.results.TransientSolution[i3].Temperature,
                md.results.TransientSolution[i3].Waterfraction,
                md.results.TransientSolution[i3].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[i3].Watercolumn,
                md.results.TransientSolution[i4].Enthalpy,
                md.results.TransientSolution[i4].Temperature,
                md.results.TransientSolution[i4].Waterfraction,
                md.results.TransientSolution[i4].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[i4].Watercolumn]
