#Test Name: SquareSheetShelfDiadSSA3dDakotaSamp
from os import getcwd
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from partitioner import *
from dmeth_params_set import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md.materials.rho_ice = 10**7  #involved in the mass flux, make it easy
md.geometry.thickness[:] = 1  #make it easy
md.geometry.surface = md.geometry.base + md.geometry.thickness

#constrain all velocities to 1 m / yr, in the y - direction
md.stressbalance.spcvx[:] = 0
md.stressbalance.spcvy[:] = 1
md.stressbalance.spcvz[:] = 0

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = 20
partition = partitioner(md, 'package', 'chaco', 'npart', npart, 'weighting', 'on')[0] - 1
md.qmu.isdakota = 1

#variables
md.qmu.variables.drag_coefficient = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_FrictionCoefficient',
    'mean', np.ones((npart, 1)),
    'stddev', .01 * np.ones((npart, 1)),
    'partition', partition)

#responses
md.qmu.responses.MaxVel = response_function.response_function('descriptor', 'MaxVel')
#md.qmu.responses.IceVolume = response_function.response_function('descriptor', 'IceVolume')  #commented for matlab consistency
md.qmu.responses.MassFlux1 = response_function.response_function('descriptor', 'indexed_MassFlux_1')
md.qmu.responses.MassFlux2 = response_function.response_function('descriptor', 'indexed_MassFlux_2')
md.qmu.responses.MassFlux3 = response_function.response_function('descriptor', 'indexed_MassFlux_3')
md.qmu.responses.MassFlux4 = response_function.response_function('descriptor', 'indexed_MassFlux_4')
md.qmu.responses.MassFlux5 = response_function.response_function('descriptor', 'indexed_MassFlux_5')
md.qmu.responses.massFlux6 = response_function.response_function('descriptor', 'indexed_MassFlux_6')
md.qmu.responses.massFlux7 = response_function.response_function('descriptor', 'indexed_MassFlux_7')

#mass flux profiles
md.qmu.mass_flux_profiles = ['../Exp/MassFlux1.exp', '../Exp/MassFlux2.exp', '../Exp/MassFlux3.exp', '../Exp/MassFlux4.exp', '../Exp/MassFlux5.exp', '../Exp/MassFlux6.exp', '../Exp/Square.exp']
md.qmu.mass_flux_profile_directory = getcwd()

# nond_sampling study
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(md.qmu.method, 'seed', 1234, 'samples', 20, 'sample_type', 'lhs')

# parameters
md.qmu.params.interval_type = 'forward'
md.qmu.params.direct = True
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

# There may be a pair of numpy warnings in the function true_divide,
#       this is normal and will not affect the results
#       See src / m / qmu / dakota_out_parse.py, function "dak_tab_out" for details
md = solve(md, 'Stressbalance', 'overwrite', 'y')

#Fields and tolerances to track changes
md.qmu.results = md.results.dakota

#ok, mass flux of 3 profiles should be-3 Gt / yr - 3 Gt / yr and the sum, which is - 6 Gt / yr
#we recover those mass fluxes through the mean of the response.
#also, we recover the max velo, which should be 1m / yr.
#we put all that data in the montecarlo field, which we will use to test for success.
#also, check that the stddev are 0.
md.results.dakota.montecarlo = []
for i in range(8):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].mean)

for i in range(8):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['moments']
field_tolerances = [1e-11]
field_values = [md.results.dakota.montecarlo]
