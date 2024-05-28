#Test Name: SquareSheetConstrainedTherTranNyeH2O
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from nye import *

md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')

md = md.extrude(3, 1.)
md = setflowequation(md, 'SSA', 'all')

#Transient options
md.cluster = generic('name', gethostname(), 'np', 3)
md.materials.rheology_law = 'NyeH2O'
md.materials.rheology_B = nye(md.initialization.temperature, 2)

md.transient.isstressbalance = 0
md.transient.ismasstransport = 0
md.transient.issmb = 1
md.transient.isthermal = 1
md.transient.isgroundingline = 0
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Temperature3', 'BasalforcingsGroundediceMeltingRate3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate]
