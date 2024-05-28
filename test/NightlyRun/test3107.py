#Test Name: SquareShelfConstrainedMasstransp3dAdolcMumps
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from issmmumpssolver import issmmumpssolver


md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.extrude(5, 3.)
md.cluster = generic('name', gethostname(), 'np', 3)

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = issmmumpssolver()
md = solve(md, 'Masstransport')

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.MasstransportSolution.Thickness]
