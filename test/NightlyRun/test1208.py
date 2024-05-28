#Test Name: EISMINTA
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


"""
EISMINT benchmark experiment A
"""

numlayers = 8
resolution = 50000.

#To begin with the numerical model
md = triangle(model(), '../Exp/SquareEISMINT750000.exp', resolution)
md = setmask(md, '', '')
md = parameterize(md, '../Par/RoundSheetEISMINT.py')

#We extrude the model to have a 3d model
md.extrude(numlayers, 1.)
md = setflowequation(md, 'SIA', 'all')

#Spc the nodes on the bed
pos = np.where(md.mesh.vertexonbase)
md.stressbalance.spcvx[pos] = 0.
md.stressbalance.spcvy[pos] = 0.
md.stressbalance.spcvz[pos] = 0.

#Adapt the time steps to the resolution
md.timestepping.time_step = 15.
md.settings.output_frequency = 500
md.timestepping.final_time = 30000.
md.masstransport.stabilization = 1
md.thermal.stabilization = 1

#Now we can solve the problem
md.cluster = generic('name', gethostname(), 'np', 8)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vz', 'Vel', 'Pressure', 'Thickness', 'Base', 'Surface', 'Temperature', 'BasalforcingsGroundediceMeltingRate']
field_tolerances = [1e-08, 1e-08, 1e-07, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-07, 1e-07]
field_values = [md.results.TransientSolution[-1].Vx,
                md.results.TransientSolution[-1].Vy,
                md.results.TransientSolution[-1].Vz,
                md.results.TransientSolution[-1].Vel,
                md.results.TransientSolution[-1].Pressure,
                md.results.TransientSolution[-1].Thickness,
                md.results.TransientSolution[-1].Base,
                md.results.TransientSolution[-1].Surface,
                md.results.TransientSolution[-1].Temperature,
                md.results.TransientSolution[-1].BasalforcingsGroundediceMeltingRate]
