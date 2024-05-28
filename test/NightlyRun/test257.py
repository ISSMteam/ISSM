#Test Name: SquareShelfSMBarma
import numpy as np
from model import *
from socket import gethostname
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from SMBarma import SMBarma
from solve import solve
from triangle import triangle


md = triangle(model(), '../Exp/Square.exp', 80000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.transient.requested_outputs = ['default', 'IceVolume', 'SmbMassBalance']

ymax = max(md.mesh.y)
xmax = max(md.mesh.x)
# Generate basin IDs for 3 basins
idbasin = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.y >= 2. / 3. * ymax)[0]
iid2 = np.where(np.logical_and(md.mesh.y < 2. / 3. * ymax, md.mesh.x >= 1. / 3. * xmax))[0]
iid3 = np.where(np.logical_and(md.mesh.y < 2. / 3. * ymax, md.mesh.x < 1. / 3. * xmax))[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idbasin[ii] = 1
    if idbasin[ii] == 0:  # no vertex was found in basin 1
        for vertex in range(3):
            if md.mesh.elements[ii][vertex] - 1 in iid2:  # one vertex in basin 2; NOTE: offset because of 1-based vertex indexing
                idbasin[ii] = 2
    if idbasin[ii] == 0:  # no vertex was found in basin 1 and 2
        idbasin[ii] = 3

# SMB parameters
numparams                  = 2
numbreaks                  = 1
intercept                  = np.array([[0.5,1.0],[1.0,0.6],[2.0,3.0]]) 
trendlin                   = np.array([[0.0,0.0],[0.01,0.001],[-0.01,0]])
polynomialparams           = np.transpose(np.stack((intercept,trendlin)),(1,2,0)) 
datebreaks                 = np.array([[3],[3],[3]]);

md.timestepping.start_time = 0
md.timestepping.time_step = 1/12
md.timestepping.final_time = 2
md.smb = SMBarma()
md.smb.num_basins = 3  # number of basins
md.smb.basin_id = idbasin  # prescribe basin ID number to elements;
md.smb.num_params        = numparams
md.smb.num_breaks        = numbreaks
md.smb.polynomialparams  = polynomialparams
md.smb.datebreaks        = datebreaks
md.smb.ar_order = 4
md.smb.ma_order = 1
md.smb.arma_timestep = 2.0  #timestep of the ARMA model [yr]
md.smb.arlag_coefs = np.array([[0.2, 0.1, 0.05, 0.01], [0.4, 0.2, -0.2, 0.1], [0.4, -0.4, 0.1, -0.1]])
md.smb.malag_coefs = np.array([[1.0],[0],[0.2]])

lm0                   = np.array([1e-4*np.array([1,-0.1,-1]),1e-6*np.array([1,-0.1,-1]),1e-5*np.array([1,-0.1,-1])])
lm1                   = np.array([1e-4*np.array([2,-0.2,-2]),1e-6*np.array([2,-0.2,-2]),1e-5*np.array([2,-0.2,-2])])
lm2                   = np.array([1e-4*np.array([3,-0.3,-3]),1e-6*np.array([3,-0.3,-3]),1e-5*np.array([3,-0.3,-3])])
lm3                   = np.array([1e-4*np.array([4,-0.4,-4]),1e-6*np.array([4,-0.4,-4]),1e-5*np.array([4,-0.4,-4])])
lm4                   = np.array([1e-4*np.array([5,-0.5,-5]),1e-6*np.array([5,-0.5,-5]),1e-5*np.array([5,-0.5,-5])])
lm5                   = np.array([1e-4*np.array([6,-0.6,-6]),1e-6*np.array([6,-0.6,-6]),1e-5*np.array([6,-0.6,-6])])
lm6                   = np.array([1e-4*np.array([7,-0.7,-7]),1e-6*np.array([7,-0.7,-7]),1e-5*np.array([7,-0.7,-7])])
lm7                   = np.array([1e-4*np.array([8,-0.8,-8]),1e-6*np.array([8,-0.8,-8]),1e-5*np.array([8,-0.8,-8])])
lm8                   = np.array([1e-4*np.array([9,-0.9,-9]),1e-6*np.array([9,-0.9,-9]),1e-5*np.array([9,-0.9,-9])])
lm9                   = np.array([1e-4*np.array([10,-1.0,-10]),1e-6*np.array([10,-1.0,-10]),1e-5*np.array([10,-1.0,-10])])
lm10                  = np.array([1e-4*np.array([11,-1.1,-11]),1e-6*np.array([11,-1.1,-11]),1e-5*np.array([11,-1.1,-11])])
lm11                  = np.array([1e-4*np.array([12,-1.2,-12]),1e-6*np.array([12,-1.2,-12]),1e-5*np.array([12,-1.2,-12])])
md.smb.lapserates     = np.stack((lm0,lm1,lm2,lm3,lm4,lm5,lm6,lm7,lm8,lm9,lm10,lm11),axis=2)
ebins                 = np.array([[100,300],[200,400],[250,450]])
md.smb.elevationbins  = np.stack([ebins for ii in range(12)],axis=2)


# Stochastic forcing
md.stochasticforcing.isstochasticforcing = 1
md.stochasticforcing.fields = ['SMBarma']
md.stochasticforcing.covariance = np.array([[0.15, 0.08, -0.02], [0.08, 0.12, -0.05], [-0.02, -0.05, 0.1]])  # global covariance among- and between-fields
md.stochasticforcing.randomflag = 0  # fixed random seeds
md.stochasticforcing.stochastictimestep  = 1.0

md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Thickness1', 'Volume1', 'SmbMassBalance1',
    'Vx2', 'Vy2', 'Vel2', 'Thickness2', 'Volume2', 'SmbMassBalance2',
    'Vx3', 'Vy3', 'Vel3', 'Thickness3', 'Volume3', 'SmbMassBalance3'
]
field_tolerances = [
    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
    1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12
]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].IceVolume,
    md.results.TransientSolution[0].SmbMassBalance,
    md.results.TransientSolution[11].Vx,
    md.results.TransientSolution[11].Vy,
    md.results.TransientSolution[11].Vel,
    md.results.TransientSolution[11].Thickness,
    md.results.TransientSolution[11].IceVolume,
    md.results.TransientSolution[11].SmbMassBalance,
    md.results.TransientSolution[23].Vx,
    md.results.TransientSolution[23].Vy,
    md.results.TransientSolution[23].Vel,
    md.results.TransientSolution[23].Thickness,
    md.results.TransientSolution[23].IceVolume,
    md.results.TransientSolution[23].SmbMassBalance
]
