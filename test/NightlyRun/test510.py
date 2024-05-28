#Test Name: PigSteaFS
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Pig.exp', 20000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md.extrude(2, 1.)
md = setflowequation(md, 'FS', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md = solve(md, 'Steadystate')

# Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-06]
field_values = [md.results.SteadystateSolution.Vx,
                md.results.SteadystateSolution.Vy,
                md.results.SteadystateSolution.Vz,
                md.results.SteadystateSolution.Vel,
                md.results.SteadystateSolution.Pressure,
                md.results.SteadystateSolution.Temperature,
                md.results.SteadystateSolution.BasalforcingsGroundediceMeltingRate]
