#Test Name: SquareShelfSMBGembClimConstrainT
import numpy as np
import sys
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from SMBgemb import *

md = triangle(model(), '../Exp/Square.exp', 350000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.materials.rho_ice = 910
md.cluster = generic('name', gethostname(), 'np', 3)

#Use of Gemb method for SMB computation
md.smb = SMBgemb(md.mesh)
md.smb.dsnowIdx = 3
md.smb.aIdx = 2
md.smb.denIdx = 1

#load hourly surface forcing date from 1979 to 2009:
if sys.version_info.major == 2:
    inputs = np.load('../Data/gemb_input.npy', allow_pickle=True).item()
else:
    inputs = np.load('../Data/gemb_input.npy', allow_pickle=True, encoding='bytes').item()

#setup the inputs:
md.smb.Ta = np.append(np.tile(np.conjugate(inputs[b'Ta0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.V = np.append(np.tile(np.conjugate(inputs[b'V0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.dswrf = np.append(np.tile(np.conjugate(inputs[b'dsw0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.dlwrf = np.append(np.tile(np.conjugate(inputs[b'dlw0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.P = np.append(np.tile(np.conjugate(inputs[b'P0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.eAir = np.append(np.tile(np.conjugate(inputs[b'eAir0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.pAir = np.append(np.tile(np.conjugate(inputs[b'pAir0']), (md.mesh.numberofelements, 1)), np.conjugate([inputs[b'dateN']]), axis=0)
md.smb.Vz = np.tile(np.conjugate(inputs[b'LP']['Vz']), (md.mesh.numberofelements, 1)).flatten()
md.smb.Tz = np.tile(np.conjugate(inputs[b'LP']['Tz']), (md.mesh.numberofelements, 1)).flatten()
md.smb.Tmean = np.tile(np.conjugate(inputs[b'LP']['Tmean']), (md.mesh.numberofelements, 1)).flatten()
md.smb.C = np.tile(np.conjugate(inputs[b'LP']['C']), (md.mesh.numberofelements, 1)).flatten()

md.smb.Ta = md.smb.Ta[:, 0:365 * 8]
md.smb.V = md.smb.V[:, 0:365 * 8]
md.smb.dswrf = md.smb.dswrf[:, 0:365 * 8]
md.smb.dlwrf = md.smb.dlwrf[:, 0:365 * 8]
md.smb.P = md.smb.P[:, 0:365 * 8]
md.smb.eAir = md.smb.eAir[:, 0:365 * 8]
md.smb.pAir = md.smb.pAir[:, 0:365 * 8]

md.timestepping.cycle_forcing = 1
md.smb.isconstrainsurfaceT = 1

#smb settings
md.smb.requested_outputs = ['SmbDz', 'SmbT', 'SmbD', 'SmbRe', 'SmbGdn', 'SmbGsp', 'SmbEC', 'SmbA', 'SmbMassBalance', 'SmbMAdd', 'SmbDzAdd', 'SmbFAC']

#only run smb core:
md.transient.isstressbalance = 0
md.transient.ismasstransport = 0
md.transient.isthermal = 0

#time stepping:
md.timestepping.start_time = 1965.6
md.timestepping.final_time = 1966.6
md.timestepping.time_step = 1. / 365.0
md.timestepping.interp_forcing = 0.

#Run transient
md = solve(md, 'Transient')

nlayers = md.results.TransientSolution[0].SmbT.shape[1]
for i in range(1, len(md.results.TransientSolution)):
    nlayers=np.minimum(md.results.TransientSolution[i].SmbT.shape[1], nlayers)

#Fields and tolerances to track changes
field_names = [
    'Layers',
    'SmbDz1', 'SmbT1', 'SmbD1', 'SmbRe1', 'SmbGdn1', 'SmbGsp1', 'SmbA1', 'SmbEC1', 'SmbMassBalance1', 'SmbMAdd1', 'SmbDzAdd1', 'SmbFAC1',
    'SmbDz2', 'SmbT2', 'SmbD2', 'SmbRe2', 'SmbGdn2', 'SmbGsp2', 'SmbA2', 'SmbEC2', 'SmbMassBalance2', 'SmbMAdd2', 'SmbDzAdd2', 'SmbFAC2',
    'SmbDz3', 'SmbT3', 'SmbD3', 'SmbRe3', 'SmbGdn3', 'SmbGsp3', 'SmbA3', 'SmbEC3', 'SmbMassBalance3', 'SmbMAdd3', 'SmbDzAdd3', 'SmbFAC3',
    'SmbDz4', 'SmbT4', 'SmbD4', 'SmbRe4', 'SmbGdn4', 'SmbGsp4', 'SmbA4', 'SmbEC4', 'SmbMassBalance4', 'SmbMAdd4', 'SmbDzAdd4', 'SmbFAC4'
]
field_tolerances = [
    1e-12,
    1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,1e-12,
    1e-12,4e-12,1e-11,1e-10,4e-11,1e-11,1e-12,1e-11,1e-12,1e-12,1e-12,1e-11,
    1e-12,4e-12,2e-12,2e-11,1e-10,1e-11,1e-12,1e-11,1e-11,1e-12,1e-12,1e-11,
    1e-11,1e-11,1e-10,1e-11,2e-11,3e-11,1e-12,4e-12,1e-10,1e-12,1e-12,2e-11
]
# Shape is different in python solution (fixed using reshape) which can cause test failure
field_values = [
    nlayers,
    md.results.TransientSolution[0].SmbDz[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbT[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbD[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbRe[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbGdn[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbGsp[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbA[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[0].SmbEC[0],
    md.results.TransientSolution[0].SmbMassBalance[0],
    md.results.TransientSolution[0].SmbMAdd[0],
    md.results.TransientSolution[0].SmbDzAdd[0],
    md.results.TransientSolution[0].SmbFAC[0],
    md.results.TransientSolution[145].SmbDz[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbT[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbD[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbRe[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbGdn[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbGsp[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbA[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[145].SmbEC[0],
    md.results.TransientSolution[145].SmbMassBalance[0],
    md.results.TransientSolution[145].SmbMAdd[0],
    md.results.TransientSolution[145].SmbDzAdd[0],
    md.results.TransientSolution[145].SmbFAC[0],
    md.results.TransientSolution[146].SmbDz[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbT[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbD[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbRe[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbGdn[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbGsp[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbA[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[146].SmbEC[0],
    md.results.TransientSolution[146].SmbMassBalance[0],
    md.results.TransientSolution[146].SmbMAdd[0],
    md.results.TransientSolution[146].SmbDzAdd[0],
    md.results.TransientSolution[146].SmbFAC[0],
    md.results.TransientSolution[-1].SmbDz[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbT[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbD[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbRe[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbGdn[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbGsp[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbA[0, 0:nlayers].reshape(1, -1),
    md.results.TransientSolution[-1].SmbEC[0],
    md.results.TransientSolution[-1].SmbMassBalance[0],
    md.results.TransientSolution[-1].SmbMAdd[0],
    md.results.TransientSolution[-1].SmbDzAdd[0],
    md.results.TransientSolution[-1].SmbFAC[0]
]
