#Test Name: ISMIPCFS
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

"""
This test is a test from the ISMP - HOM Intercomparison project.
Pattyn and Payne 2006
"""

L_list = [80000]
results = []

for L in L_list:
    md = triangle(model(), "../Exp/Square_{}.exp".format(L), L / 10.)  #size 3 * L
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPC.py')
    md.friction.coefficient = np.sqrt(md.constants.yts * (1000. + 1000. * np.sin(md.mesh.x * 2. * np.pi / L) * np.sin(md.mesh.y * 2. * np.pi / L)))
    md.extrude(10, 1.)

#Add spc on the borders
    pos = np.where(np.logical_or.reduce((md.mesh.x == 0., md.mesh.x == np.max(md.mesh.x), md.mesh.y == 0., md.mesh.y == np.max(md.mesh.y))))
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.
    if (L == 5000.):
        md.stressbalance.spcvx[pos] = 15.66
        md.stressbalance.spcvy[pos] = -0.1967
    elif (L == 10000.):
        md.stressbalance.spcvx[pos] = 16.04
        md.stressbalance.spcvy[pos] = -0.1977
    elif (L == 20000.):
        md.stressbalance.spcvx[pos] = 16.53
        md.stressbalance.spcvy[pos] = -1.27
    elif (L == 40000.):
        md.stressbalance.spcvx[pos] = 17.23
        md.stressbalance.spcvy[pos] = -3.17
    elif (L == 80000.):
        md.stressbalance.spcvx[pos] = 16.68
        md.stressbalance.spcvy[pos] = -2.69
    elif (L == 160000.):
        md.stressbalance.spcvx[pos] = 16.03
        md.stressbalance.spcvy[pos] = -1.27

    md = setflowequation(md, 'FS', 'all')

#Compute the stressbalance
    md.cluster = generic('name', gethostname(), 'np', 8)
    md = solve(md, 'Stressbalance')

#Plot the results and save them
    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    results.append(md.results.StressbalanceSolution)

#   plotmodel(md, 'data', vx, 'data', vy, 'data', vz, 'layer  #all', md.mesh.numberoflayers)

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
