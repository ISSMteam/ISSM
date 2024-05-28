#Test Name: PigTranRignotarma 
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
md.timestepping.time_step = 0.05
md.timestepping.final_time = 2

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
#Dimensionalities
nb_tf = 2

# Calving parameters
md.mask.ice_levelset = 1e4 * (md.mask.ice_levelset + 0.5)
md.calving.calvingrate = 0 * np.ones((md.mesh.numberofvertices,))
md.levelset.spclevelset = np.full((md.mesh.numberofvertices,), np.nan)
md.levelset.migration_max = 10.0
### Frontal forcing parameters ###
md.frontalforcings = frontalforcingsrignotarma()
# Polynomial params #
numparams        = 3;
numbreaks        = 2;
intercept        = np.array([[2.5,2.0,0.1],[0.5,0.5,1.5]])
trendlin         = np.array([[-0.5,-0.2,0.1],[0,0,0]])
trendquad        = np.array([[0,0.0,0],[0.1,0.1,0.1]])
datebreaks       = np.array([[0.5,1.5],[0.5,1.5]])
polynomialparams = np.transpose(np.stack((intercept,trendlin,trendquad)),(1,2,0))
# Monthly effects params #
numbreaksM       = 1;
intcpsMp0        = np.array([[-0.5,-0.5,0,0,0,0,0.5,0.5,0,0,0,0],
                    [-1.0,-1.0,0,0,0,0,1.0,1.0,0,0,0,0]])
intcpsMp1        = np.array([[-0.25,-0.25,0,0,0,0,0.,0.,0,0,0,0],
                    [-0.1,-0.1,0,0,0,0,0.1,0.1,0,0,0,0]])
intcpsM          = np.transpose(np.stack((intcpsMp0,intcpsMp1)),(1,2,0)) 
trendsMp0        = np.array([[0,0,0,0,0,0,0.,0.0,0,0,0,0],
                    [0.0,0.0,0,-0.0,0,0,0.0,0.0,0,0,0,0]])
trendsMp1        = np.array([[0,-0.12,0,0,0,0,0.,0.0,0,0.0,0,0],
                    [0.0,-0.1,0,-0.0,0,0,0.0,0.0,0,0,0,0]])
trendsM          = np.transpose(np.stack((trendsMp0,trendsMp1)),(1,2,0)) 
datebreaksM      = np.array([[1],[1]]) 
# Subglacial discharge params #
isdischargearma            = 1
sd_ar_order                = 1
sd_ma_order                = 1
sd_num_breaks              = 1
sd_num_params              = 2
sd_arma_timestep           = 1
sd_arlag_coefs             = np.array([[0.95],[0.95]])
sd_malag_coefs             = np.array([[0.0],[0.0]])
sd_datebreaks              = np.array([[1.0],[1.0]])
sd_monthlyfrac             = np.array([[0,0,0,0,0,0,0.5,0.5,0,0,0,0],[0,0,0,0,0,0,0.5,0.5,0,0,0,0]])
sd_const                   = np.array([[50000,70000],[8000,10000.0]])
sd_trend                   = np.array([[0.0,10000],[0,0]])
sd_polyparam               = np.transpose(np.stack((sd_const,sd_trend)),(1,2,0))


md.frontalforcings.num_basins = nb_tf
md.frontalforcings.basin_id = idb_tf
md.frontalforcings.subglacial_discharge = 0.01 * np.ones((md.mesh.numberofvertices,))
md.frontalforcings.num_params = numparams #number of parameters in the polynomial
md.frontalforcings.num_breaks = numbreaks #number of breakpoints
md.frontalforcings.polynomialparams = polynomialparams
md.frontalforcings.datebreaks = datebreaks
md.frontalforcings.ar_order = 4
md.frontalforcings.ma_order = 2
md.frontalforcings.arma_timestep = 2  # timestep of the ARMA model [yr]
md.frontalforcings.arlag_coefs = np.array([[0.1, -0.1, 0.01, -0.01], [0.2, -0.2, 0.1, 0.0]])  # autoregressive parameters
md.frontalforcings.malag_coefs = np.array([[0.1, 0.0], [0.0, 0.1]])  # moving-average parameters
md.frontalforcings.monthlyvals_numbreaks   = numbreaksM
md.frontalforcings.monthlyvals_intercepts  = intcpsM
md.frontalforcings.monthlyvals_trends      = trendsM
md.frontalforcings.monthlyvals_datebreaks  = datebreaksM
md.frontalforcings.isdischargearma         = isdischargearma
if(isdischargearma==0):
   md.frontalforcings.subglacial_discharge    = 0.01*ones(md.mesh.numberofvertices,1)
else:
    md.frontalforcings.sd_num_breaks         = sd_num_breaks
    md.frontalforcings.sd_num_params         = sd_num_params
    md.frontalforcings.sd_ar_order           = sd_ar_order
    md.frontalforcings.sd_ma_order           = sd_ma_order
    md.frontalforcings.sd_arma_timestep      = sd_arma_timestep
    md.frontalforcings.sd_arlag_coefs        = sd_arlag_coefs
    md.frontalforcings.sd_malag_coefs        = sd_malag_coefs
    md.frontalforcings.sd_datebreaks         = sd_datebreaks
    md.frontalforcings.sd_monthlyfrac        = sd_monthlyfrac
    md.frontalforcings.sd_polynomialparams   = sd_polyparam
#Floating Ice Melt parameters
md.basalforcings.floatingice_melting_rate = 0.0 * np.ones((md.mesh.numberofvertices,))

#Covariance matrix
covtf = 1e-4 * np.identity(nb_tf)
covsd = 1e3 * np.identity(nb_tf)
#covglob          = np.zeros([6,6])
#covglob[0:2,0:2] = covtf
#covglob[2:4,2:4] = covclv
#covglob[4:6,4:6] = covflmlt

#Hard-coding covariance matrix because python is complaining
covglob = np.array([[1e-4, 0., 0., 0.],
                    [0., 1e-4, 0., 0.],
                    [0., 0., 1e3, 0.],
                    [0., 0., 0., 1e3]])
#testchol = np.linalg.cholesky(covglob)
#print(testchol)

# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['FrontalForcingsRignotarma','FrontalForcingsSubglacialDischargearma']
md.stochasticforcing.defaultdimension = 2
md.stochasticforcing.default_id = idb_tf
md.stochasticforcing.covariance = covglob  # global covariance among- and between-fields
md.stochasticforcing.randomflag = 0  # determines true/false randomness

md.transient.ismovingfront = 1
md.transient.isgroundingline = 1

md.transient.requested_outputs = ['default', 'CalvingCalvingrate', 'CalvingMeltingrate', 'BasalforcingsFloatingiceMeltingRate']
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Thickness1', 'MaskIceLevelset1', 'CalvingMeltingrate1',
    'Vx2', 'Vy2', 'Vel2', 'Thickness2', 'MaskIceLevelset2', 'CalvingMeltingrate2',
    'Vx10', 'Vy10', 'Vel10', 'Thickness10', 'MaskIceLevelset10', 'CalvingMeltingrate10']

field_tolerances = [
        1e-11,2e-11,2e-11,1e-11,1e-9,1e-10,
        2e-11,1e-11,1e-11,9e-11,2e-9,1e-10,
        2e-6,1e-6,1e-6,1e-6,5e-6,1e-6]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].MaskIceLevelset,
    md.results.TransientSolution[0].CalvingMeltingrate,
    md.results.TransientSolution[19].Vx,
    md.results.TransientSolution[19].Vy,
    md.results.TransientSolution[19].Vel,
    md.results.TransientSolution[19].Thickness,
    md.results.TransientSolution[19].MaskIceLevelset,
    md.results.TransientSolution[19].CalvingMeltingrate,
    md.results.TransientSolution[39].Vx,
    md.results.TransientSolution[39].Vy,
    md.results.TransientSolution[39].Vel,
    md.results.TransientSolution[39].Thickness,
    md.results.TransientSolution[39].MaskIceLevelset,
    md.results.TransientSolution[39].CalvingMeltingrate]
