#Test Name: SquareSheetShelfTranSSA3d
import numpy as np
from model import *
from socket import gethostname

from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from frictioncoulomb import frictioncoulomb

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.geometry.bed = copy.deepcopy(md.geometry.base)
pos = np.nonzero(md.mask.ocean_levelset < 0.)
md.geometry.bed[pos] = md.geometry.bed[pos] - 10
md.friction = frictioncoulomb()
md.friction.coefficient = 20 * np.ones(md.mesh.numberofvertices)
md.friction.p = 1 * np.ones(md.mesh.numberofelements)
md.friction.q = 1 * np.ones(md.mesh.numberofelements)
md.friction.coefficientcoulomb = 0.02 * np.ones(md.mesh.numberofvertices)
md.transient.isthermal = False
md.transient.isgroundingline = True
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.transient.requested_outputs = ['default', 'GroundedArea', 'FloatingArea', 'TotalFloatingBmb', 'TotalGroundedBmb', 'TotalSmb']
md.extrude(3, 1.)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = [
    'Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1', 'TotalSmb1', 
    'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'GroundedArea2', 'FloatingArea2', 'TotalFloatingBmb2', 'TotalGroundedBmb2', 'TotalSmb2', 
    'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'GroundedArea3', 'FloatingArea3', 'TotalFloatingBmb3', 'TotalGroundedBmb3', 'TotalSmb3'
]
field_tolerances = [
    2e-13, 2e-13, 2e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
    2e-13, 2e-13, 2e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
    2e-13, 2e-13, 2e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 
]
field_values = [
    md.results.TransientSolution[0].Vx,
    md.results.TransientSolution[0].Vy,
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Pressure,
    md.results.TransientSolution[0].Base,
    md.results.TransientSolution[0].Surface,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].GroundedArea,
    md.results.TransientSolution[0].FloatingArea,
    md.results.TransientSolution[0].TotalFloatingBmb,
    md.results.TransientSolution[0].TotalGroundedBmb,
    md.results.TransientSolution[0].TotalSmb,
    md.results.TransientSolution[1].Vx,
    md.results.TransientSolution[1].Vy,
    md.results.TransientSolution[1].Vel,
    md.results.TransientSolution[1].Pressure,
    md.results.TransientSolution[1].Base,
    md.results.TransientSolution[1].Surface,
    md.results.TransientSolution[1].Thickness,
    md.results.TransientSolution[1].GroundedArea,
    md.results.TransientSolution[1].FloatingArea,
    md.results.TransientSolution[1].TotalFloatingBmb,
    md.results.TransientSolution[1].TotalGroundedBmb,
    md.results.TransientSolution[1].TotalSmb,
    md.results.TransientSolution[2].Vx,
    md.results.TransientSolution[2].Vy,
    md.results.TransientSolution[2].Vel,
    md.results.TransientSolution[2].Pressure,
    md.results.TransientSolution[2].Base,
    md.results.TransientSolution[2].Surface,
    md.results.TransientSolution[2].Thickness,
    md.results.TransientSolution[2].GroundedArea,
    md.results.TransientSolution[2].FloatingArea,
    md.results.TransientSolution[2].TotalFloatingBmb,
    md.results.TransientSolution[2].TotalGroundedBmb,
    md.results.TransientSolution[2].TotalSmb
]
