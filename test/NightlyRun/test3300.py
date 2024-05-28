import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from transient import *
from setflowequation import *
from solve import *


from generic import generic

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md.transient = transient.deactivateall(md.transient)
md.transient.ishydrology = True
md.transient.issmb = True
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md.hydrology = hydrologydc()
md.hydrology = md.hydrology.initialize(md)

md.hydrology.isefficientlayer = 1
md.hydrology.sedimentlimit_flag = 1
md.hydrology.sedimentlimit = 400.0
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))
md.hydrology.sediment_thickness = 20.0
md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = 2.0 * np.ones((md.mesh.numberofvertices))
md.hydrology.sediment_transmitivity = 1.5e-4 * np.ones((md.mesh.numberofvertices))

md.initialization.epl_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_thickness = np.ones((md.mesh.numberofvertices))
md.hydrology.spcepl_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.mask_eplactive_node = np.zeros((md.mesh.numberofvertices))

md.hydrology.epl_conductivity = 1.5e-2
md.hydrology.epl_initial_thickness = 1.0
md.hydrology.epl_colapse_thickness = 1.0e-6
md.hydrology.epl_thick_comp = 1
md.hydrology.epl_max_thickness = 5.0

md.hydrology.transfer_flag = 1.0
md.hydrology.leakage_factor = 3.9e-13

times = np.arange(0, 8.001, 0.002)
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices + 1, len(times)))

md.basalforcings.groundedice_melting_rate[:, np.where(times <= 6.0)] = -0.2
md.basalforcings.groundedice_melting_rate[:, np.where(times <= 1.0)] = 1.0
md.basalforcings.groundedice_melting_rate[-1, :] = times

md.timestepping.time_step = 0.002
md.timestepping.final_time = 8.0

md = solve(md, 'Transient')

# sedvol = np.zeros(4000)
# eplvol = np.zeros(4000)
# totvol = np.zeros(4001)
# time = np.arange(0.002, 8.001, 0.002)
# store = md.constants.g * md.hydrology.sediment_porosity * md.materials.rho_freshwater * ((md.hydrology.sediment_compressibility / md.hydrology.sediment_porosity) + md.hydrology.water_compressibility)
# sedstore = 20.0 * store
# for i in range(0, 4000):
#       sedvol[i]=np.mean(md.results.TransientSolution[i].SedimentHead)*sedstore
#       eplvol[i]=np.mean(md.results.TransientSolution[i].EplHead)*store*np.mean(md.results.TransientSolution[i].HydrologydcEplThicknessSubstep)
#       totvol[i+1]=totvol[i]+md.basalforcings.groundedice_melting_rate[0, i]*0.002

field_names = ['SedimentWaterHead5', 'EplWaterHead5', 'SedimentWaterHead40', 'EplWaterHead40']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[5].SedimentHead,
                md.results.TransientSolution[5].EplHead,
                md.results.TransientSolution[40].SedimentHead,
                md.results.TransientSolution[40].EplHead]
