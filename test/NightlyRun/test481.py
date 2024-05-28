#Test Name: TransientFrictionSchoof
import numpy as np

from frictionschoof import frictionschoof
from socket import gethostname
from model import *
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from solve import solve
from transient import transient
from triangle import triangle

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.extrude(4, 1)
md = setflowequation(md, 'HO', 'all')
md.transient.isthermal = 0
md.friction = frictionschoof(md.friction)
md.friction.C    = pow(20.e4, 0.5) * np.ones((md.mesh.numberofvertices, 1))
md.friction.Cmax = 0.5 * np.ones((md.mesh.numberofvertices, 1))
md.friction.m    = 1./3.* np.ones((md.mesh.numberofelements, 1))
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1',
    'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2',
    'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3'
]
field_tolerances = [
    2e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09,
    1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09,
    2e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09, 1e-09
]
field_values = [
    md.results.TransientSolution[0].Vx,
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
    md.results.TransientSolution[2].Thickness
]
