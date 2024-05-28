#Test Name: SquareNoDynExtrudedHydrologyDCTwoLayers
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from solve import solve
from generic import generic

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareNoDyn.py')
md.cluster = generic('name', gethostname(), 'np', 1)

md.transient.ishydrology = True
md.hydrology = hydrologydc()
md.hydrology = md.hydrology.initialize(md)

md.hydrology.isefficientlayer = 1
md.hydrology.sedimentlimit_flag = 1
md.hydrology.transfer_flag = 0
md.hydrology.sedimentlimit = 800.0
md.hydrology.mask_thawed_node = np.ones((md.mesh.numberofvertices))
md.initialization.sediment_head = np.zeros((md.mesh.numberofvertices))
md.hydrology.spcsediment_head = np.nan * np.ones((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = 2.0 * np.ones((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.hydrology.sediment_transmitivity = 3 * np.ones((md.mesh.numberofvertices))

md.initialization.epl_head = np.zeros((md.mesh.numberofvertices))
md.initialization.epl_thickness = np.ones((md.mesh.numberofvertices))
md.hydrology.spcepl_head = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.mask_eplactive_node = np.zeros((md.mesh.numberofvertices))
md.hydrology.epl_conductivity = 30
md.hydrology.epl_initial_thickness = 1
md.hydrology.epl_colapse_thickness = 1.0e-3
md.hydrology.epl_thick_comp = 1
md.hydrology.epl_max_thickness = 1
md.timestepping.time_step = 0.2
md.timestepping.final_time = 2.0

md.extrude(3, 1.)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['SedimentWaterHead1', 'EplWaterHead1', 'SedimentHeadResidual1',
               'SedimentWaterHead4', 'EplWaterHead4', 'SedimentHeadResidual4',
               'SedimentWaterHead5', 'EplWaterHead5', 'SedimentHeadResidual5',
               'SedimentWaterHead9', 'EplWaterHead9', 'SedimentHeadResidual9']
field_tolerances = [1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 5e-12, 2e-11,
                    1e-13, 5e-12, 2e-11]
field_values = [md.results.TransientSolution[0].SedimentHead,
                md.results.TransientSolution[0].EplHead,
                md.results.TransientSolution[0].SedimentHeadResidual,
                md.results.TransientSolution[3].SedimentHead,
                md.results.TransientSolution[3].EplHead,
                md.results.TransientSolution[3].SedimentHeadResidual,
                md.results.TransientSolution[4].SedimentHead,
                md.results.TransientSolution[4].EplHead,
                md.results.TransientSolution[4].SedimentHeadResidual,
                md.results.TransientSolution[8].SedimentHead,
                md.results.TransientSolution[8].EplHead,
                md.results.TransientSolution[8].SedimentHeadResidual]
