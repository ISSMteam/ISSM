#Test Name: SquareShelfDamageEvolutionSSA2dPralong
import numpy as np
from triangle import triangle
from model import *
from socket import gethostname
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from generic import generic
from solve import solve
from matdamageice import matdamageice

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md.materials = matdamageice()
md = parameterize(md, '../Par/SquareShelf.py')
md.damage.isdamage = 1
md.damage.D = 0.1 * np.ones(md.mesh.numberofvertices)
md.damage.spcdamage = np.nan * np.ones(md.mesh.numberofvertices)
md.damage.law = 1

md.damage.c1 = 1.e-11
md.damage.c2 = 0.4
md.damage.c3 = 1.e-3
md.damage.healing = 0.4
md.damage.stress_threshold = 1.e5
md.damage.stabilization = 1

md.damage.requested_outputs = ['default', 'DamageF']

md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'DamageEvolution')

field_names = ['D', 'F']
field_tolerances = [1.e-13, 1.e-13]
field_values = [md.results.DamageEvolutionSolution.DamageDbar,
                md.results.DamageEvolutionSolution.DamageF]
