#Test Name: PigTranCalvingDevdHO3d
import numpy as np
from calvingvonmises import *
from socket import gethostname
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *

md = triangle(model(), '../Exp/Pig.exp', 10000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md.extrude(5, 1.)
md = setflowequation(md, 'HO', 'all')
md.timestepping.time_step = 2
md.timestepping.final_time = 50

#calving parameters
md.mask.ice_levelset = 1e4 * (md.mask.ice_levelset + 0.5)
md.calving = calvingvonmises()
md.frontalforcings.meltingrate = np.zeros((md.mesh.numberofvertices, ))
md.transient.ismovingfront = 1
md.levelset.spclevelset = np.nan * np.ones((md.mesh.numberofvertices, ))
pos = np.where(md.mesh.vertexonboundary)
md.levelset.spclevelset[pos] = md.mask.ice_levelset[pos]
md.levelset.migration_max = 1e10
md.transient.requested_outputs = ['default', 'IceVolume', 'IceVolumeAboveFloatation','TotalSmb','TotalGroundedBmb','TotalFloatingBmb']

#Force MUMPS sequential analysis
md.toolkits.DefaultAnalysis.mat_mumps_icntl_28 = 1
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'MaskIceLevelset1',
               'IceVolume1', 'IceVolumeAboveFloatation1', 'TotalSmb1', 'TotalGroundedBmb1', 'TotalFloatingBmb1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'MaskIceLevelset2',
               'IceVolume2', 'IceVolumeAboveFloatation2', 'TotalSmb2', 'TotalGroundedBmb2', 'TotalFloatingBmb2',
               'Vx10', 'Vy10', 'Vel10', 'Pressure10', 'Bed10', 'Surface10', 'Thickness10', 'MaskIceLevelset10',
               'IceVolume10', 'IceVolumeAboveFloatation10', 'TotalSmb10', 'TotalGroundedBmb10', 'TotalFloatingBmb10']
field_tolerances = [1e-11, 2e-11, 2e-11, 1e-12, 2e-11, 6e-12, 9e-12, 2e-12,
                    1e-11, 2e-11, 2e-11, 3e-9, 2e-11,
                    2e-11, 1e-11, 1e-11, 9e-12, 2e-11, 3e-11, 2e-11, 1e-11,
                    1e-11, 2e-11, 2e-11, 8e-08, 2e-11,
                    2e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-9,
                    1e-11, 2e-11, 2e-11, 8e-08, 2e-11]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskIceLevelset,
                md.results.TransientSolution[0].IceVolume,
                md.results.TransientSolution[0].IceVolumeAboveFloatation,
                md.results.TransientSolution[0].TotalSmb,
                md.results.TransientSolution[0].TotalGroundedBmb,
                md.results.TransientSolution[0].TotalFloatingBmb,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskIceLevelset,
                md.results.TransientSolution[1].IceVolume,
                md.results.TransientSolution[1].IceVolumeAboveFloatation,
                md.results.TransientSolution[1].TotalSmb,
                md.results.TransientSolution[1].TotalGroundedBmb,
                md.results.TransientSolution[1].TotalFloatingBmb,
                md.results.TransientSolution[9].Vx,
                md.results.TransientSolution[9].Vy,
                md.results.TransientSolution[9].Vel,
                md.results.TransientSolution[9].Pressure,
                md.results.TransientSolution[9].Base,
                md.results.TransientSolution[9].Surface,
                md.results.TransientSolution[9].Thickness,
                md.results.TransientSolution[9].MaskIceLevelset,
                md.results.TransientSolution[9].IceVolume,
                md.results.TransientSolution[9].IceVolumeAboveFloatation,
                md.results.TransientSolution[9].TotalSmb,
                md.results.TransientSolution[9].TotalGroundedBmb,
                md.results.TransientSolution[9].TotalFloatingBmb]
