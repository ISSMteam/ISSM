#Test Name: EISMINTStress2
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

#test 5 and 6:
md = model()
md = triangle(md, '../Exp/SquareEISMINT.exp', 5100.)  #test3
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

vx = md.results.StressbalanceSolution.Vx
vy = md.results.StressbalanceSolution.Vy

#plot results
#plotmodel(md, 'data', vx, 'contourlevels', {0, 20, 40, 60, 80, 100, -20, -40, -60, -80, -100}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag2vx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag2vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])
#plotmodel(md, 'data', vy, 'contourlevels', { -100, -200, -300, -400, -500, -600, -700, -800, -900, -1000}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag2vy', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag2vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy']
field_tolerances = [1e-13, 1e-13]
field_values = [vx, vy]
