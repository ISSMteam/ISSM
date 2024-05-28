#Test Name: SquareShelfConstrainedStressSSA2dAdolcMumps
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from issmmumpssolver import issmmumpssolver


md = triangle(model(), '../Exp/Square.exp', 50000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.stressbalance.requested_outputs = ['default', 'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = issmmumpssolver()
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel', 'Pressure',
               'DeviatoricStressxx', 'DeviatoricStressyy', 'DeviatoricStressxy']
field_tolerances = [1e-12, 1e-12, 1e-12, 1e-12,
                    1e-12, 1e-12, 1e-12]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.DeviatoricStressxx,
                md.results.StressbalanceSolution.DeviatoricStressyy,
                md.results.StressbalanceSolution.DeviatoricStressxy]
