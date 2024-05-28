#Test Name: SquareShelfTranForceNeg2dDakotaSampLinearPart
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
from dmeth_params_set import *

md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

md.timestepping.time_step = 1
md.settings.output_frequency = 1
md.timestepping.final_time = 4

smb = np.ones((md.mesh.numberofvertices, )) * 3.6
smb = np.array([smb, smb * -1]).T

md.smb.mass_balance = smb
md.smb.mass_balance = np.concatenate((md.smb.mass_balance, [[1.5, 3]]))
md.transient.isthermal = 0
#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = md.mesh.numberofvertices
partition = partitioner(md, 'package', 'linear', 'npart', npart)[0] - 1

#variables
md.qmu.variables.surface_mass_balance = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_SmbMassBalance',
    'mean', np.ones((npart, 1)),
    'stddev', .1 * np.ones((npart, 1)),
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

#mass flux profiles
md.qmu.mass_flux_profiles = ['../Exp/MassFlux1.exp', '../Exp/MassFlux2.exp', '../Exp/MassFlux3.exp', '../Exp/MassFlux4.exp', '../Exp/MassFlux5.exp', '../Exp/MassFlux6.exp']
md.qmu.mass_flux_profile_directory = getcwd()

#  nond_sampling study
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(md.qmu.method, 'seed', 1234, 'samples', 20, 'sample_type', 'lhs')
dver = str(version)
if ((int(dver[0]) == 4 and int(dver[2]) > 2) or int(dver[0]) > 4):
    md.qmu.method = dmeth_params_set(md.qmu.method, 'rng', 'rnum2')

#parameters
md.qmu.params.direct = True
md.qmu.params.analysis_components = ''
md.qmu.params.interval_type = 'forward'
md.qmu.params.tabular_graphics_data = True
md.qmu.isdakota = 1

if version >= 6:
    md.qmu.params.analysis_driver = 'matlab'
    md.qmu.params.evaluation_scheduling = 'master'
    md.qmu.params.processors_per_evaluation = 2
else:
    md.qmu.params.analysis_driver = 'stressbalance'
    md.qmu.params.evaluation_concurrency = 1

md.stressbalance.reltol = 10**-5  #tighten for qmu analyses
md.transient.requested_outputs = ['IceVolume']

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Transient', 'overwrite', 'y')
md.qmu.results = md.results.dakota

#Fields and tolerances to track changes
md.results.dakota.moments = []
for i in range(8):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].mean)

for i in range(8):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['moments']
field_tolerances = [1e-11]
field_values = [md.results.dakota.moments]
