#Test Name: SquareSheetShelfDiadSSA3dDakotaMassFlux
import numpy as np
from os import getcwd
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from partitioner import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.materials.rho_ice = 10**7  #involved in the mass flux, make it easy
md.geometry.thickness[:] = 1  #make it easy
md.geometry.surface = md.geometry.base + md.geometry.thickness

#constrain all velocities to 1 m / yr, in the y - direction
md.stressbalance.spcvx = np.zeros((md.mesh.numberofvertices, ))
md.stressbalance.spcvy = np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvz = np.zeros((md.mesh.numberofvertices, ))

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = 20
partition = partitioner(md, 'package', 'chaco', 'npart', npart, 'weighting', 'on')[0] - 1

#variables
md.qmu.variables.drag_coefficient = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_FrictionCoefficient',
    'mean', np.ones((npart, 1)),
    'stddev', .01 * np.ones((npart, 1)),
    'partition', partition
    )

#responses
md.qmu.responses.MaxVel = response_function.response_function('descriptor','MaxVel')
md.qmu.responses.IceVolume = response_function.response_function('descriptor','IceVolume')
md.qmu.responses.MassFlux1 = response_function.response_function('descriptor','indexed_MassFlux_1')
md.qmu.responses.MassFlux2 = response_function.response_function('descriptor','indexed_MassFlux_2')
md.qmu.responses.MassFlux3 = response_function.response_function('descriptor','indexed_MassFlux_3')
md.qmu.responses.MassFlux4 = response_function.response_function('descriptor','indexed_MassFlux_4')
md.qmu.responses.MassFlux5 = response_function.response_function('descriptor','indexed_MassFlux_5')
md.qmu.responses.massFlux6 = response_function.response_function('descriptor','indexed_MassFlux_6')
md.qmu.responses.massFlux7 = response_function.response_function('descriptor','indexed_MassFlux_7')

#mass flux profiles
md.qmu.mass_flux_profiles = ['../Exp/MassFlux1.exp', '../Exp/MassFlux2.exp', '../Exp/MassFlux3.exp', '../Exp/MassFlux4.exp', '../Exp/MassFlux5.exp', '../Exp/MassFlux6.exp', '../Exp/Square.exp']
md.qmu.mass_flux_profile_directory = getcwd()

#method
md.qmu.method = dakota_method.dakota_method('nond_l')

#parameters
md.qmu.params.direct = True
md.qmu.params.interval_type = 'forward'
md.qmu.isdakota = 1
md.stressbalance.reltol = 10**-5  #tighten for qmu analyses

if version >= 6:
    md.qmu.params.analysis_driver = 'matlab'
    md.qmu.params.evaluation_scheduling = 'master'
    md.qmu.params.processors_per_evaluation = 2
else:
    md.qmu.params.analysis_driver = 'stressbalance'
    md.qmu.params.evaluation_concurrency = 1

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Stressbalance', 'overwrite', 'y')
md.qmu.results = md.results.dakota

#Fields and tolerances to track changes
#ok, mass flux of 3 profiles should be-3 Gt / yr - 3 Gt / yr and the sum, which is - 6 Gt / yr
#we recover those mass fluxes through the mean of the response.
#also, we recover the max velo, which should be 1m / yr.
#we put all that data in the moments, which we will use to test for success.
#also, check that the stddev are 0.
md.results.dakota.moments = []
for i in range(8):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].mean)

for i in range(8):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['moments']
field_tolerances = [1e-11]
field_values = [md.results.dakota.moments]
