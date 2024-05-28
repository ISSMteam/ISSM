#Test Name: PigTherTranSUPG
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Pig.exp', 30000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')
md.thermal.stabilization = 2
md.cluster = generic('name', gethostname(), 'np', 3)
md.transient.isstressbalance = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = False
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Temperature1', 'BasalforcingsGroundediceMeltingRate1',
               'Temperature2', 'BasalforcingsGroundediceMeltingRate2']
field_tolerances = [1e-13, 1e-8, 1e-13, 5e-8]
field_values = [md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate]
