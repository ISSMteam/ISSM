#Test Name: EISMINTStress1
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

#Compute solution for SSA's model
md.cluster = generic('name', gethostname(), 'np', 8)
md = solve(md, 'Stressbalance')

#plot results
vx = md.results.StressbalanceSolution.Vx
vy = md.results.StressbalanceSolution.Vy

#plotmodel(md, 'data', vx, 'contourlevels', {0, 20, 40, 60, 60, 100, 120, 140, 160, 180, -20, -40, -60, -80, -100, -120, -140, -160, -180}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag1vx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag1vx.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#plotmodel(md, 'data', vy, 'contourlevels', { -100, -200, -300, -400, -500, -600, -700, -800, -900, -1000}, ...
#       'contourcolor', 'k')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('eismintdiag1vy', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 2, 'hardcopy', 'off')
#       system(['mv eismintdiag1vy.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceShelf '])

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy']
field_tolerances = [1e-13, 1e-13]
field_values = [vx, vy]
