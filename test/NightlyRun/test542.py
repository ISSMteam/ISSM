#Test Name: PigTranFrontalforcingsrignot
import numpy as np
from calvingvonmises import calvingvonmises
from frontalforcingsrignot import frontalforcingsrignot
from socket import gethostname
from model import *
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from solve import solve
from triangle import triangle


md = triangle(model(), '../Exp/Pig.exp', 10000)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md = setflowequation(md, 'SSA', 'all')
md.timestepping.time_step = 2
md.timestepping.final_time = 50

# Separate domain in 2 basins
idbasin = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.x <= -1.6e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idbasin[ii] = 1
    if idbasin[ii] == 0:  # no vertex was found in basin 1
        idbasin[ii] = 2

# Calving and frontalforcings parameters
md.mask.ice_levelset = 1e4 * (md.mask.ice_levelset + 0.5)
md.calving = calvingvonmises()
md.frontalforcings = frontalforcingsrignot()
md.frontalforcings.num_basins = 2
md.frontalforcings.basin_id = idbasin
md.frontalforcings.subglacial_discharge = 0.1 * np.ones((md.mesh.numberofvertices,))
md.frontalforcings.thermalforcing = 0.5 * np.ones((md.mesh.numberofvertices,))
for elem in range(md.mesh.numberofelements):
    if idbasin[elem] == 2:
        md.frontalforcings.thermalforcing[md.mesh.elements[elem, 0:3] - 1] = 1.5  #NOTE: offset because of 1-based vertex indexing

md.transient.ismovingfront = 1
md.levelset.spclevelset = np.full((md.mesh.numberofvertices,), np.nan)
md.levelset.migration_max = 1e10

md.transient.requested_outputs = ['default', 'CalvingMeltingrate']
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'MaskIceLevelset1', 'CalvingMetlingRate1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'MaskIceLevelset2', 'CalvingMetlingRate2',
               'Vx10', 'Vy10', 'Vel10', 'Pressure10', 'Bed10', 'Surface10', 'Thickness10', 'MaskIceLevelset10', 'CalvingMetlingRate10']
field_tolerances = [
    1e-11, 2e-11, 2e-11, 1e-12, 2e-11, 6e-12, 9e-12, 1e-12, 1e-9,
    2e-11, 1e-11, 1e-11, 9e-12, 2e-1, 2e-11, 2e-11, 1e-11, 1e-9,
    2e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-9, 1e-9]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Pressure,
    md.results.TransientSolution[0].Base,
    md.results.TransientSolution[0].Surface,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].MaskIceLevelset,
    md.results.TransientSolution[0].CalvingMeltingrate,
    md.results.TransientSolution[1].Vx,
    md.results.TransientSolution[1].Vy,
    md.results.TransientSolution[1].Vel,
    md.results.TransientSolution[1].Pressure,
    md.results.TransientSolution[1].Base,
    md.results.TransientSolution[1].Surface,
    md.results.TransientSolution[1].Thickness,
    md.results.TransientSolution[1].MaskIceLevelset,
    md.results.TransientSolution[1].CalvingMeltingrate,
    md.results.TransientSolution[9].Vx,
    md.results.TransientSolution[9].Vy,
    md.results.TransientSolution[9].Vel,
    md.results.TransientSolution[9].Pressure,
    md.results.TransientSolution[9].Base,
    md.results.TransientSolution[9].Surface,
    md.results.TransientSolution[9].Thickness,
    md.results.TransientSolution[9].MaskIceLevelset,
    md.results.TransientSolution[9].CalvingMeltingrate
]
