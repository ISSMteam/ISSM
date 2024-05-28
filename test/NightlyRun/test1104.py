#Test Name: ISMIPBFS
import numpy as np
from model import *
from socket import gethostname
from squaremesh import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

"""
This test is a test from the ISMP - HOM Intercomparison project.
Pattyn and Payne 2006
"""

L_list = [80000.]
results = []

for L in L_list:
    nx = 20  #numberof nodes in x direction
    ny = 20
    md = model()
    md = squaremesh(md, L, L, nx, ny)
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPB.py')
    md.extrude(10, 1.)
    md = setflowequation(md, 'HO', 'all')

#Create dirichlet on the bed only
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

    pos = np.where(md.mesh.vertexonbase)
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.

#Create MPCs to have periodic boundary conditions
    posx = np.where(md.mesh.x == 0.)[0]
    posx2 = np.where(md.mesh.x == np.max(md.mesh.x))[0]

    posy = np.where(np.logical_and.reduce((md.mesh.y == 0., md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]  #Don't take the same nodes two times
    posy2 = np.where(np.logical_and.reduce((md.mesh.y == np.max(md.mesh.y), md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]

    md.stressbalance.vertex_pairing = np.vstack((np.vstack((posx + 1, posx2 + 1)).T, np.vstack((posy + 1, posy2 + 1)).T))

#Compute the stressbalance
    md.stressbalance.abstol = np.nan
    md.cluster = generic('name', gethostname(), 'np', 8)
    md = solve(md, 'Stressbalance')
    pos = np.where(np.logical_or.reduce((md.mesh.x == 0., md.mesh.y == 0., md.mesh.x == np.max(md.mesh.x), md.mesh.y == np.max(md.mesh.y))))
    md.stressbalance.spcvx[pos] = md.results.StressbalanceSolution.Vx[pos][:, 0]
    md.stressbalance.spcvy[pos] = md.results.StressbalanceSolution.Vy[pos][:, 0]
    md.stressbalance.vertex_pairing = np.empty((0, 2), int)
    md = setflowequation(md, 'FS', 'all')
    md = solve(md, 'Stressbalance')

#Plot the results and save them
    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    results.append(md.results.StressbalanceSolution)

#       plotmodel(md, 'data', vx, 'data', vy, 'data', vz, 'layer#all', md.mesh.numberoflayers)

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-08, 1e-07, 1e-08]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
