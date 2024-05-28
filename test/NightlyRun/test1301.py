#Test Name: ThermalMelting
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from paterson import *
from solve import *


"""
This file can be run to check that the melting in simple conduction is correctly modeled.
There is no velocity (no advection) the only thermal boundary conditions are an imposed temperature
at upper surface and an impose flux at its base. The result must be a linear temperature from the upper to the lower
surface with an imposed slope (Geothermal flux). if it is not the case, something is thermal modeling has been changed...
"""

printingflag = False

md = model()
md = triangle(md, '../Exp/Square.exp', 100000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareThermal.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')

#Some conditions specific to melting test
md.initialization.pressure = np.zeros((md.mesh.numberofvertices, 1), int)
md.initialization.temperature = 273.15 * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(md.mesh.vertexonsurface)[0]
md.thermal.spctemperature[pos] = md.initialization.temperature[pos]
md.materials.rheology_B = paterson(md.initialization.temperature)

#analytical results
#melting heat = geothermal flux
#Mb * L * rho = G = > Mb = G / L * rho
melting = md.basalforcings.geothermalflux / (md.materials.rho_ice * md.materials.latentheat) * md.constants.yts

#modeled results
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Thermal')

#plot results
comp_melting = md.results.ThermalSolution.BasalforcingsGroundediceMeltingRate
relative = np.abs((comp_melting - melting) / melting) * 100.
relative[np.nonzero(comp_melting == melting)[0]] = 0.
#plotmodel(md, 'data', comp_melting, 'title', 'Modeled melting', 'data', melting, 'title', 'Analytical melting', ...
#       'data', comp_melting - melting, 'title', 'Absolute error', 'data', relative, 'title', 'Relative error [%]', ...
#       'layer  #all', 1, 'caxis  #2', [1.02964 1.02966] * 1.0e - 4, 'FontSize  #all', 20, 'figposition', 'mathieu')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('thermalmelting', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 0.7, 'hardcopy', 'off')
#       system(['mv thermalmelting.png ' ISSM_DIR '/website/doc_pdf/validation/Images/Thermal '])

#Fields and tolerances to track changes
field_names = ['BasalMelting']
field_tolerances = [1e-08]
field_values = [comp_melting]
