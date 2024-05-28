#Test Name: SquareSheetShelfSteaEnthalpyHO3dDakotaSampNeff

# TODO:
# - Figure out why test fails intermittently on Mac with "IndexError: list 
# index out of range"
#

from os import getcwd
from socket import gethostname

import numpy as np

from ContourToMesh import *
from dmeth_params_set import *
from model import *
from parameterize import *
from partitioner import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *


#model not consistent:  equality thickness = surface-base violated

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.extrude(3, 2.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.timestepping.time_step = 0.
md.thermal.isenthalpy = 1
md.thermal.isdynamicbasalspc = 1
md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices, 1))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices, 1))

md.friction.coupling = 3
md.friction.effective_pressure = md.materials.rho_ice * md.constants.g * md.geometry.thickness + md.materials.rho_water * md.constants.g * md.geometry.base

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = 10
partition, md = partitioner(md, 'package', 'chaco', 'npart', npart, 'weighting', 'on')
partition -= 1
md.qmu.isdakota = 1

#variables
md.qmu.variables.neff = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_FrictionEffectivePressure',
    'mean', np.ones((npart, 1)),
    'stddev', .05 * np.ones((npart, 1)),
    'partition', partition
    )
md.qmu.variables.geoflux = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_BasalforcingsGeothermalflux',
    'mean', np.ones((npart, 1)),
    'stddev', .05 * np.ones((npart, 1)),
    'partition', partition
    )

#responses
md.qmu.responses.MaxVel = response_function.response_function('descriptor','MaxVel')
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
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(
    md.qmu.method,
    'seed', 1234,
    'samples', 20,
    'sample_type', 'random'
    )

#parameters
md.qmu.params.direct = True
md.qmu.params.analysis_components = ''
md.qmu.params.tabular_graphics_data = True

if version >= 6:
    md.qmu.params.analysis_driver = 'matlab'
    md.qmu.params.evaluation_scheduling = 'master'
    md.qmu.params.processors_per_evaluation = 2
else:
    md.qmu.params.analysis_driver = 'stressbalance'
    md.qmu.params.evaluation_concurrency = 1

md.stressbalance.reltol = 10**-5  #tighten for qmu analyses

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Steadystate', 'overwrite', 'y')

#Fields and tolerances to track changes
md.qmu.results = md.results.dakota

#we put all the mean and stdev data in the montecarlo field, which we will use to test for success.
md.results.dakota.montecarlo = []
for i in range(8):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].mean)

for i in range(8):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['montecarlo']
field_tolerances = [2e-10]
field_values = [md.results.dakota.montecarlo]
