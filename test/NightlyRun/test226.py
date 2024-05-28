#Test Name: SquareShelfTranCflSSA2d
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping = timesteppingadaptive(md.timestepping)
md.timestepping.final_time = 15.
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1',
               'Bed1', 'Surface1', 'Thickness1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2',
               'Bed2', 'Surface2', 'Thickness2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3',
               'Bed3', 'Surface3', 'Thickness3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness]
