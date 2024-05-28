#Test Name: 79NorthStochFrictionWaterPressure
import numpy as np

from socket import gethostname
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *


md = triangle(model(), '../Exp/79North.exp', 6000)
md = setmask(md, '../Exp/79NorthShelf.exp', '')
md = parameterize(md, '../Par/79North.py')
md = setflowequation(md, 'SSA', 'all')

#Default friction
md.friction = friction()
md.friction.coefficient = 30 * np.ones(md.mesh.numberofvertices)
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Basin separation default
idb_df = np.zeros((md.mesh.numberofelements))
iid1 = np.where(md.mesh.y <= -1.08e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idb_df[ii] = 1
    if idb_df[ii] == 0:  # no vertex was found in basin 1
        for vertex in range(3):
            idb_df[ii] = 2
#Covariance matrix
covPw = np.array([[0.75e10, 0.0], [0.0, 0.5e10]])

# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['FrictionWaterPressure']
md.stochasticforcing.defaultdimension = 2
md.stochasticforcing.default_id = idb_df
md.stochasticforcing.covariance = covPw  # global covariance
md.stochasticforcing.randomflag = 0  # determines true/false randomness

md.transient.issmb = 0
md.transient.ismasstransport = 1
md.transient.isstressbalance = 1
md.transient.isthermal = 0
md.transient.isgroundingline = 0

md.transient.requested_outputs = ['default', 'FrictionWaterPressure']
md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 5
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Thickness1', 'FrictionWaterPressure1',
               'Vx2', 'Vy2', 'Vel2', 'Thickness2', 'FrictionWaterPressure2',
               'Vx10', 'Vy10', 'Vel10', 'Thickness10', 'FrictionWaterPressure10']

field_tolerances = [2e-10, 2e-10, 2e-10, 2e-10, 2e-10,
                    4e-10, 4e-10, 4e-10, 4e-10, 4e-10,
                    8e-10, 8e-10, 8e-10, 8e-10, 8e-10]

field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].FrictionWaterPressure,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].FrictionWaterPressure,
                md.results.TransientSolution[4].Vx,
                md.results.TransientSolution[4].Vy,
                md.results.TransientSolution[4].Vel,
                md.results.TransientSolution[4].Thickness,
                md.results.TransientSolution[4].FrictionWaterPressure]
