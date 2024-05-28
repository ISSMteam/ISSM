#Test Name: PigTranStochasticforcingCovariance
import numpy as np
from frontalforcingsrignotarma import frontalforcingsrignotarma
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
md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 10

# Basin separation TF
idb_tf = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.x <= -1.6e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idb_tf[ii] = 1
    if idb_tf[ii] == 0:  # no vertex was found in basin 1
        for vertex in range(3):
            idb_tf[ii] = 2
# Basin separation default
idb_df = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.x <= -1.62e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idb_df[ii] = 1
    if idb_df[ii] == 0:  # no vertex was found in basin 1
        for vertex in range(3):
            idb_df[ii] = 2
#Dimensionalities
nb_tf = 2
nb_clv  = 2
nb_flmlt = 2

# Calving parameters
md.mask.ice_levelset = 1e4 * (md.mask.ice_levelset + 0.5)
md.calving.calvingrate = 0.3 * np.ones((md.mesh.numberofvertices,))
md.levelset.spclevelset = np.full((md.mesh.numberofvertices,), np.nan)
md.levelset.migration_max = 10.0
### Frontal forcing parameters ###
md.frontalforcings = frontalforcingsrignotarma()
md.frontalforcings.num_basins = nb_tf
md.frontalforcings.basin_id = idb_tf
# Polynomial params #
md.frontalforcings.num_params       = 1 #only a constant term
md.frontalforcings.num_breaks       = 0 #no breakpoint
constval                            = np.array([[2.5],[0.5]])
md.frontalforcings.polynomialparams = np.copy(constval)
# No monthly effects: do nothing #
# ARMA model parameters #
md.frontalforcings.ar_order = 3
md.frontalforcings.ma_order = 2
md.frontalforcings.arma_timestep = 2  # timestep of the ARMA model [yr]
md.frontalforcings.arlag_coefs = np.array([[0.1, -0.1, 0.01], [0.2, -0.2, 0.1]])  # autoregressive parameters
md.frontalforcings.malag_coefs = np.array([[0.1, 0.0], [0.0, 0.1]])  # moving-average parameters
# No ARMA model of subglacial discharge: simply specify values at vertices #
md.frontalforcings.subglacial_discharge = 10 * np.ones((md.mesh.numberofvertices,))

#Floating Ice Melt parameters
md.basalforcings.floatingice_melting_rate = 0.1 * np.ones((md.mesh.numberofvertices,))

#Covariance matrix
covtf = 1e-4 * np.identity(nb_tf)
covclv = 1e-1 * np.identity(nb_clv)
covclv[0, 0] = 1 / 10 * covclv[0, 0]
covflmlt = 0.05 * np.identity(nb_flmlt)
#covglob          = np.zeros([6,6])
#covglob[0:2,0:2] = covtf
#covglob[2:4,2:4] = covclv
#covglob[4:6,4:6] = covflmlt

#Hard-coding covariance matrix because python is complaining
covglob = np.array([[1e-4, 0., 0., 0., 0., 0.],
                    [0., 1e-4, 0., 0., 0., 0.],
                    [0., 0., 1e-2, 0., 0., 0.],
                    [0., 0., 0., 1e-1, 0., 0.],
                    [0., 0., 0., 0., 0.05, 0.],
                    [0., 0., 0., 0., 0., 0.05]])
#testchol = np.linalg.cholesky(covglob)
#print(testchol)

# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['FrontalForcingsRignotarma', 'DefaultCalving', 'FloatingMeltRate']
md.stochasticforcing.defaultdimension = 2
md.stochasticforcing.default_id = idb_df
md.stochasticforcing.covariance = covglob  # global covariance among- and between-fields
md.stochasticforcing.randomflag = 0  # determines true/false randomness

md.transient.ismovingfront = 1
md.transient.isgroundingline = 1

md.transient.requested_outputs = ['default', 'CalvingCalvingrate', 'CalvingMeltingrate', 'BasalforcingsFloatingiceMeltingRate']
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Thickness1', 'MaskIceLevelset1', 'CalvingCalvingrate1', 'CalvingMeltingrate1', 'BasalforcingsFloatingiceMeltingRate1',
    'Vx5', 'Vy5', 'Vel5', 'Thickness5', 'MaskIceLevelset5', 'CalvingCalvingrate5', 'CalvingMeltingrate5', 'BasalforcingsFloatingiceMeltingRate5',
    'Vx10', 'Vy10', 'Vel10', 'Thickness10', 'MaskIceLevelset10', 'CalvingCalvingrate10', 'CalvingMeltingrate10', 'BasalforcingsFloatingiceMeltingRate10']

field_tolerances = [
    1e-11, 2e-11, 2e-11, 1e-11, 1e-9, 1e-10, 1e-10, 1e-10,
    2e-11, 1e-11, 1e-11, 9e-11, 2e-9, 1e-10, 1e-10, 1e-10,
    2e-6, 1e-6, 1e-6, 1e-6, 5e-6, 1e-6, 1e-6, 1e-6]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].MaskIceLevelset,
    md.results.TransientSolution[0].CalvingCalvingrate,
    md.results.TransientSolution[0].CalvingMeltingrate,
    md.results.TransientSolution[0].BasalforcingsFloatingiceMeltingRate,
    md.results.TransientSolution[4].Vx,
    md.results.TransientSolution[4].Vy,
    md.results.TransientSolution[4].Vel,
    md.results.TransientSolution[4].Thickness,
    md.results.TransientSolution[4].MaskIceLevelset,
    md.results.TransientSolution[4].CalvingCalvingrate,
    md.results.TransientSolution[4].CalvingMeltingrate,
    md.results.TransientSolution[4].BasalforcingsFloatingiceMeltingRate,
    md.results.TransientSolution[9].Vx,
    md.results.TransientSolution[9].Vy,
    md.results.TransientSolution[9].Vel,
    md.results.TransientSolution[9].Thickness,
    md.results.TransientSolution[9].MaskIceLevelset,
    md.results.TransientSolution[9].CalvingCalvingrate,
    md.results.TransientSolution[9].CalvingMeltingrate,
    md.results.TransientSolution[9].BasalforcingsFloatingiceMeltingRate]
