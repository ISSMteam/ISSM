#Test Name: SquareShelfConstrainedTranAdolcReverseVsForward
import copy
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from independent import *
from dependent import *
from SetIceShelfBC import *
from solve import *


#This test runs test3020 with autodiff on, and checks that
#the value of the scalar forward difference match a step - wise differential

#First configure
md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md.transient.requested_outputs = ['IceVolume', 'MaxVel']
md.verbose = verbose('autodiff', True)
md.stressbalance.restol = 0.000001

#setup autodiff parameters
index = 1  #this is the scalar component we are checking against
md.autodiff.independents = [independent('name', 'md.geometry.thickness', 'type', 'vertex', 'nods', md.mesh.numberofvertices, 'fos_forward_index', index)]
md.autodiff.dependents = [dependent('name', 'IceVolume', 'type', 'scalar'),
                          dependent('name', 'MaxVel', 'type', 'scalar')]
md.autodiff.driver = 'fos_forward'

#PYTHON: indices start at 0, make sure to offset index
index = index - 1

#parameters for the step - wise derivative
delta = 0.00001
h1 = md.geometry.thickness[index]
h0 = h1 * (1. - delta)
h2 = h1 * (1. + delta)
deltaH = (h2 - h0)

#save model:
md2 = copy.deepcopy(md)

#evaluate derivative by forward and backward stepping
#forward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h0
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = SetIceShelfBC(md)

md = solve(md, 'Transient')
V0 = md.results.TransientSolution[-1].IceVolume
MaxV0 = md.results.TransientSolution[-1].MaxVel

#backward
md = copy.deepcopy(md2)
md.autodiff.isautodiff = False
md.geometry.thickness[index] = h2
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = SetIceShelfBC(md)

md = solve(md, 'Transient')
V2 = md.results.TransientSolution[-1].IceVolume
MaxV2 = md.results.TransientSolution[-1].MaxVel

#compute resulting derivative
dVdh_an = (V2 - V0) / deltaH
dMaxVdh_an = (MaxV2 - MaxV0) / deltaH

#evaluate derivative using ADOLC
md = copy.deepcopy(md2)
md.autodiff.isautodiff = True
md.geometry.thickness[index] = h1
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = SetIceShelfBC(md)

md = solve(md, 'Transient')
#retrieve directly
dVdh_ad = md.results.TransientSolution[0].AutodiffJacobian[0]
dMaxVdh_ad = md.results.TransientSolution[0].AutodiffJacobian[1]

print("dV / dh: analytical:  %16.16g\n       using adolc:  %16.16g\n" % (dVdh_an, dVdh_ad))
print("dMaxV / dh: analytical:  %16.16g\n       using adolc:  %16.16g\n" % (dMaxVdh_an, dMaxVdh_ad))

#Fields and tolerances to track changes
field_names = ['dV/dh - dV / dh0', 'dMaxV/dh - dMaxV / dh0']
field_tolerances = [1e-13, 1e-13]
field_values = [dVdh_ad - dVdh_an, dMaxVdh_an - dMaxVdh_ad]
