#Test Name: ISMIPE
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from squaremesh import *

#This test is a test from the ISMP - HOM Intercomparison project.
#TestE
#Four tests to run: - Pattyn frozen
# - Stokes frozen
# - Pattyn with some sliding
# - Stokes with some sliding
printingflag = False
results = []

for i in range(4):
    Lx = 10.  #in m
    Ly = 5000.  #in m
    nx = 3  #number of nodes in x direction
    ny = 51
    md = model()
    md = squaremesh(md, Lx, Ly, nx, ny)
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPE.py')
    md = md.extrude(10, 1.)

    if i == 0 or i == 2:
        md = setflowequation(md, 'HO', 'all')
    elif i == 1 or i == 3:
        md = setflowequation(md, 'FS', 'all')

#Create MPCs to have periodic boundary conditions
    posx = np.where(md.mesh.x == 0.)[0]
    posx2 = np.where(md.mesh.x == max(md.mesh.x))[0]
    md.stressbalance.vertex_pairing = np.column_stack((posx, posx2))

#Create spcs on the bed
    pos = np.where(md.mesh.vertexonbase)[0]
    md.stressbalance.spcvx = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    md.stressbalance.spcvy = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    md.stressbalance.spcvz = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.
    md.stressbalance.spcvz[pos] = 0.

#Remove the spc where there is some sliding (case 3 and 4):
    if i == 2 or i == 3:
        pos = np.intersect1d(np.where((md.mesh.y / max(md.mesh.y)) >= 0.44), np.where((md.mesh.y / max(md.mesh.y)) <= 0.5))[0]
        md.stressbalance.spcvx[pos] = float('NaN')
        md.stressbalance.spcvy[pos] = float('NaN')
        md.stressbalance.spcvz[pos] = float('NaN')

#Compute the stressbalance
    md.cluster = generic('name', gethostname(), 'np', 8)
    md = solve(md, 'Stressbalance')

    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    results[i] = md.results.StressbalanceSolution


#Fields and tolerances to track changes
field_names = ['VyPattynSliding', 'VzPattynSliding',
               'VxStokesSliding', 'VyStokesSliding', 'VzStokesSliding',
               'VyPattynFrozen', 'VzPattynFrozen',
               'VxStokesFrozen', 'VyStokesFrozen', 'VzStokesFrozen']
field_tolerances = [1e-05, 1e-05,
                    1e-05, 1e-06, 1e-06,
                    1e-05, 1e-04,
                    1e-05, 1e-05, 1e-06]

field_values = []
for i in range(4):
    result = results[i]
    field_values += [result.Vx, result.Vy, result.Vz]
