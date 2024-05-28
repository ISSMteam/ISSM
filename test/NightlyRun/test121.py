#Test Name: SquareShelfConstrainedEnthalpyTran
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 180000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md.extrude(3, 1.)
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.transient.isstressbalance = False
md.transient.ismasstransport = False
md.transient.issmb = True
md.transient.isthermal = True
md.transient.isgroundingline = False
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Enthalpy1', 'Waterfraction1', 'Temperature1',
               'Enthalpy2', 'Waterfraction2', 'Temperature2',
               'Enthalpy3', 'Waterfraction3', 'Temperature3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-9, 1e-13]
field_values = [md.results.TransientSolution[0].Enthalpy,
                md.results.TransientSolution[0].Waterfraction,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[1].Enthalpy,
                md.results.TransientSolution[1].Waterfraction,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[2].Enthalpy,
                md.results.TransientSolution[2].Waterfraction,
                md.results.TransientSolution[2].Temperature]
