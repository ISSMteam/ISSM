#Test Name: EISMINTTran2
import numpy as np
import sys
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


"""
Test on the stressbalance model and the masstransport in 2d
"""

printingflag = False

#tests 3 and 4: using Glen's flow law
md = model()
md = triangle(md, '../Exp/SquareEISMINT.exp', 3550.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareEISMINT.py')
md = setflowequation(md, 'SSA', 'all')  #SSA's model and 2d

#Impose a non zero velocity on the upper boundary condition (y = max(y))
pos = np.where(md.mesh.y == np.max(md.mesh.y))
heavyside = np.where(np.logical_and(md.mesh.y == np.max(md.mesh.y), ((1. + sys.float_info.epsilon) * np.ones((np.size(md.mesh.y))) - ((md.mesh.x - 100000.) / 25000.)**2) > 0))
md.stressbalance.spcvy[pos] = np.zeros((np.size(pos)))
md.stressbalance.spcvy[heavyside] = 400. * (((md.mesh.x[heavyside] - 100000.) / 25000.)**2 - np.ones((np.size(heavyside))))

#Compute solution for SSA's model
md.cluster = generic('name', gethostname(), 'np', 8)
md = solve(md, 'Stressbalance')

#plot results
md.initialization.vx = md.results.StressbalanceSolution.Vx
md.initialization.vy = md.results.StressbalanceSolution.Vy

md.timestepping.time_step = 1.
md.timestepping.final_time = 5000.
md.masstransport.stabilization = 1
md = solve(md, 'Transient')

#plotmodel(md, 'data', (md.results.TransientSolution(end).Vx))
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eisminttrans2vx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eisminttrans2vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#plotmodel(md, 'data', (md.results.TransientSolution(end).Vy))
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eisminttrans2vy', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eisminttrans2vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#plotmodel(md, 'data', (md.results.TransientSolution(end).Thickness))
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eisminttrans2thickness', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eisminttrans2thickness.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Thickness']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[-1].Vx,
                md.results.TransientSolution[-1].Vy,
                md.results.TransientSolution[-1].Thickness]
