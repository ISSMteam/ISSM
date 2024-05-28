#Test Name: ThermalAdvection
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


"""
This file can be run to check that the advection - diffusion  is correctly modeled.
There is u = v = 0 and w = cst everywhere the only thermal boundary conditions are an imposed temperature
at upper surface and an impose flux at its base.
"""

printingflag = False

md = model()
md = triangle(md, '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareThermal.py')
md.extrude(30, 1.)  #NB: the more one extrudes, the better (10 -> relative~0.35%, 20 -> 0.1%, 30 -> 0.05%)
md = setflowequation(md, 'HO', 'all')

#Thermal boundary conditions
pos1 = np.where(np.isnan(md.mesh.lowerelements))[0]
md.thermal.spctemperature[md.mesh.elements[pos1, 0:3] - 1] = 10.
pos2 = np.where(np.isnan(md.mesh.upperelements))[0]
md.thermal.spctemperature[md.mesh.elements[pos2, 3:6] - 1] = 0.
md.initialization.vz = 0.1 * np.ones((md.mesh.numberofvertices))
md.initialization.vel = np.sqrt(md.initialization.vx**2 + md.initialization.vy**2 + md.initialization.vz**2)
md.initialization.pressure = np.zeros((md.mesh.numberofvertices), int)

md.thermal.stabilization = 2
#analytical results
#d2T / dz2 - w * rho_ice * c / k * dT / dz = 0   T(surface)=0  T(bed)=10 = > T = A exp(alpha z) + B
alpha = 0.1 / md.constants.yts * md.materials.rho_ice * md.materials.heatcapacity / md.materials.thermalconductivity  #alpha = w rho_ice c / k  and w = 0.1m / an
A = 10. / (np.exp(alpha * (-1000.)) - 1.)  #A = T(bed) / (exp(alpha * bed) - 1)  with bed= -1000 T(bed)=10
B = -A
md.initialization.temperature = A * np.exp(alpha * md.mesh.z) + B

#modeled results
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Thermal')

#plot results
comp_temp = md.results.ThermalSolution.Temperature
relative = np.abs((comp_temp - md.initialization.temperature) / md.initialization.temperature) * 100.
relative[np.nonzero(comp_temp == md.initialization.temperature)[0]] = 0.
#plotmodel(md, 'data', comp_temp, 'title', 'Modeled temperature [K]', 'data', md.initialization.temperature, 'view', 3, ...
#       'title', 'Analytical temperature [K]', 'view', 3, 'data', comp_temp - md.initialization.temperature, ...
#       'title', 'Absolute error [K]', 'view', 3, 'data', relative, 'title', 'Relative error [%]', 'view', 3, ...
#       'figposition', 'mathieu', 'FontSize  #all', 20)
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('thermaladvection', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 0.7, 'hardcopy', 'off')
#       system(['mv thermaladvection.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT '])

#Fields and tolerances to track changes
field_names = ['AdvectionTemperature']
field_tolerances = [1e-13]
field_values = [comp_temp]
