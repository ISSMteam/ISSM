#Test Name: SquareShelfConstrainedStressSSA2dAdolc
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from issmgslsolver import issmgslsolver


md = triangle(model(), '../Exp/Square.exp', 50000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']
md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = issmgslsolver()
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure',
               'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.DeviatoricStressxx,
                md.results.StressbalanceSolution.DeviatoricStressyy,
                md.results.StressbalanceSolution.DeviatoricStressxy]
