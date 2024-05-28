#Test Name: SquareShelfConstrainedMasstransp2dDGAdolcMumps
from model import *
from socket import gethostname
from triangle import *
from meshconvert import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from issmmumpssolver import issmmumpssolver


md = triangle(model(), '../Exp/Square.exp', 150000.)
md = meshconvert(md)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.masstransport.stabilization = 3
md.masstransport.spcthickness = md.geometry.thickness

md.autodiff.isautodiff = True
md.toolkits.DefaultAnalysis = issmmumpssolver()
md = solve(md, 'Masstransport')

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-13]
field_values = [md.results.MasstransportSolution.Thickness]
