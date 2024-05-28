#Test Name: SquareSheetConstrainedStressMOLHO2d
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from SetMOLHOBC import SetMOLHOBC


md = triangle(model(), '../Exp/Square.exp', 200000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md = setflowequation(md, 'MOLHO', 'all')
md = SetMOLHOBC(md)
md.extrude(5, 1.)
md.cluster = generic('name', gethostname(), 'np', 3)
md.stressbalance.requested_outputs = ['default', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']
md = solve(md, 'Stressbalance')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'VxSurface', 'VySurface', 'VxShear', 'VyShear', 'VxBase', 'VyBase']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.StressbalanceSolution.Vx,
                md.results.StressbalanceSolution.Vy,
                md.results.StressbalanceSolution.Vz,
                md.results.StressbalanceSolution.Vel,
                md.results.StressbalanceSolution.Pressure,
                md.results.StressbalanceSolution.VxSurface,
                md.results.StressbalanceSolution.VySurface,
                md.results.StressbalanceSolution.VxShear,
                md.results.StressbalanceSolution.VyShear,
                md.results.StressbalanceSolution.VxBase,
                md.results.StressbalanceSolution.VyBase]
