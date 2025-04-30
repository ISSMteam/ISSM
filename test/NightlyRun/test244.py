#Test Name: SquareShelfSMBGembDakota
import numpy as np
import scipy.io as spio
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from SMBgemb import *
from IssmConfig import *

from partitioner import *
from dakota_method import *
from normal_uncertain import *
from uniform_uncertain import *
from response_function import *
from dmeth_params_set import *

md = triangle(model(), '../Exp/Square.exp', 300000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.materials.rho_ice = 910
md.cluster = generic('name', gethostname(), 'np', 3)
md.geometry.bed = md.geometry.base

# Use of Gemb method for SMB computation
md.smb = SMBgemb(md.mesh)
md.smb.dsnowIdx = 0
md.smb.swIdx = 1

#load hourly surface forcing date from 1979 to 2009:
inputs = spio.loadmat('../Data/gemb_input.mat', squeeze_me=True)

#setup the inputs:
md.smb.Ta = np.append(np.tile(np.conjugate(inputs['Ta0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.V = np.append(np.tile(np.conjugate(inputs['V0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.dswrf = np.append(np.tile(np.conjugate(inputs['dsw0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.dlwrf = np.append(np.tile(np.conjugate(inputs['dlw0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.P = np.append(np.tile(np.conjugate(inputs['P0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.eAir = np.append(np.tile(np.conjugate(inputs['eAir0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.pAir = np.append(np.tile(np.conjugate(inputs['pAir0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.pAir = np.append(np.tile(np.conjugate(inputs['pAir0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs['dateN']]), axis=0)
md.smb.Vz = np.tile(np.conjugate(inputs['LP']['Vz']), (md.mesh.numberofelements, 1)).flatten()
md.smb.Tz = np.tile(np.conjugate(inputs['LP']['Tz']), (md.mesh.numberofelements, 1)).flatten()
md.smb.Tmean = np.tile(np.conjugate(inputs['LP']['Tmean']), (md.mesh.numberofelements, 1)).flatten()
md.smb.C = np.tile(np.conjugate(inputs['LP']['C']), (md.mesh.numberofelements, 1)).flatten()

#smb settings
md.smb.requested_outputs = ['SmbDz', 'SmbT', 'SmbD', 'SmbRe', 'SmbGdn', 'SmbGsp', 'SmbEC', 'SmbA', 'SmbMassBalance']

#only run smb core:
md.transient.isstressbalance = 0
md.transient.ismasstransport = 1
md.transient.isthermal = 0

#time stepping:
md.timestepping.start_time = 1965.
md.timestepping.final_time = 1965.75
md.timestepping.time_step = 1. / 365.0
md.timestepping.interp_forcing = 0.

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = md.mesh.numberofelements
partition = partitioner(md, 'package', 'linear', 'type', 'element', 'npart', npart)[0] - 1

#variables
md.qmu.variables.surface_mass_balance = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_SmbC',
    'mean', np.ones((npart, 1)),
    'stddev', .5 * np.ones((npart, 1)),
    'partition', partition
    )
Tmin = 273.
telms = np.atleast_2d(np.min(md.smb.Ta[0:-1, :], 1))
mint_on_partition = telms.flatten()
for pa in range(np.size(mint_on_partition)):
    vi = np.where(partition == pa)
    mint = telms[0, vi] * 1.05
    pos = np.where(mint < Tmin)
    mint[pos] = Tmin
    mint_on_partition[pa] = max(mint / telms[0, vi])

mint_on_partition[np.where(np.isnan(mint_on_partition))] = 10**-10
upper = np.maximum(np.minimum(np.maximum(1.05, mint_on_partition), 0.9999), 0.0001)
upper = upper.reshape(npart, 1)
md.qmu.variables.surface_mass_balanceTa = uniform_uncertain.uniform_uncertain(
    'descriptor', 'scaled_SmbTa',
    'lower', .95 * np.ones((npart, 1)),
    'upper', upper,
    'partition', partition
)

#responses
md.qmu.responses.IceVolume = response_function.response_function('descriptor','IceVolume')
md.qmu.responses.IceMass = response_function.response_function('descriptor','IceMass')
md.qmu.responses.TotalSmb = response_function.response_function('descriptor','TotalSmb')

#  nond_sampling study
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(md.qmu.method, 'seed', 1234, 'samples', 3, 'sample_type', 'lhs')
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
md.transient.requested_outputs = ['IceVolume', 'TotalSmb', 'IceMass']

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Transient', 'overwrite', 'y')
md.qmu.results = md.results.dakota

#Fields and tolerances to track changes
md.results.dakota.moments = []
for i in range(3):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].mean)

for i in range(3):
    md.results.dakota.moments.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['moments']
field_tolerances = [2e-6]
field_values = [md.results.dakota.moments]
