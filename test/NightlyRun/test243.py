#Test Name: SquareShelfSMBGemb
from socket import gethostname
import sys
import numpy as np
from model import *
from parameterize import *
from setflowequation import *
from setmask import *
from SMBgemb import *
from solve import *
from triangle import *

md = triangle(model(), '../Exp/Square.exp', 350000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.materials.rho_ice = 910
md.cluster = generic('name', gethostname(), 'np', 3)

#Use of Gemb method for SMB computation
md.smb = SMBgemb(md.mesh, md.geometry)
md.smb.dsnowIdx = 1
md.smb.swIdx = 1

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

#smb settings
md.smb.requested_outputs = ['SmbDz','SmbT','SmbD','SmbRe','SmbGdn','SmbGsp','SmbEC',
                            'SmbA','SmbMassBalance','SmbMAdd','SmbDzAdd','SmbFAC','SmbMeanSHF','SmbMeanLHF',
                            'SmbMeanULW','SmbNetLW','SmbNetSW','SmbWAdd','SmbRunoff','SmbRefreeze','SmbMelt',
                            'SmbEC','SmbPrecipitation','SmbRain','SmbAccumulatedMassBalance','SmbAccumulatedRunoff',
                            'SmbAccumulatedMelt','SmbAccumulatedEC','SmbAccumulatedPrecipitation','SmbAccumulatedRain',
                            'SmbAccumulatedPrecipitation','SmbAccumulatedRefreeze']

#only run smb core:
md.transient.isstressbalance = 0
md.transient.ismasstransport = 0
md.transient.isthermal = 0

#time stepping:
md.timestepping.start_time = 1965.
md.timestepping.final_time = 1966.
md.timestepping.time_step = 1.0 / 365
md.timestepping.interp_forcing = 0.

#Run transient
md = solve(md, 'Transient')

nlayers = md.results.TransientSolution[0].SmbT.shape[1]
for i in range(1, len(md.results.TransientSolution)):
    nlayers=np.minimum(md.results.TransientSolution[i].SmbT.shape[1], nlayers)

#Fields and tolerances to track changes
field_names = ['Layers', 'SmbDz', 'SmbT', 'SmbD', 'SmbRe', 'SmbGdn', 'SmbGsp', 'SmbA', 'SmbEC', 'SmbMassBalance', 'SmbMAdd', 'SmbDzAdd', 'SmbFAC', 'SmbMeanSHF', 'SmbMeanLHF', 'SmbMeanULW', 'SmbNetLW', 'SmbNetSW', 'SmbAccumulatedMassBalance', 'SmbAccumulatedRunoff', 'SmbAccumulatedMelt', 'SmbAccumulatedEC', 'SmbAccumulatedPrecipitation', 'SmbAccumulatedRain', 'SmbAccumulatedRefreeze', 'SmbRunoff', 'SmbMelt', 'SmbEC', 'SmbPrecipitation', 'SmbRain', 'SmbRefreeze', 'SmbWAdd']
field_tolerances = [1e-12, 4e-11, 2e-11, 3e-11, 6e-11, 8e-11, 8e-11, 1e-12, 5e-11, 2e-12, 1e-12, 1e-12, 4e-11, 2e-11, 5e-11, 1e-11, 9e-10, 2e-11, 1e-11, 9e-10, 2e-11, 2e-09, 1e-11, 1e-11, 1e-11, 8e-10, 2e-11, 2e-11, 1e-11, 1e-11, 1e-11, 1e-11]

# Shape is different in python solution (fixed using reshape) which can cause test failure
field_values = [
    nlayers,
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
    md.results.TransientSolution[-1].SmbFAC[0],
    md.results.TransientSolution[-1].SmbMeanSHF[0],
    md.results.TransientSolution[-1].SmbMeanLHF[0],
    md.results.TransientSolution[-1].SmbMeanULW[0],
    md.results.TransientSolution[-1].SmbNetLW[0],
    md.results.TransientSolution[-1].SmbNetSW[0],
    md.results.TransientSolution[-1].SmbAccumulatedMassBalance[0],
    md.results.TransientSolution[-1].SmbAccumulatedRunoff[0],
    md.results.TransientSolution[-1].SmbAccumulatedMelt[0],
    md.results.TransientSolution[-1].SmbAccumulatedEC[0],
    md.results.TransientSolution[-1].SmbAccumulatedPrecipitation[0],
    md.results.TransientSolution[-1].SmbAccumulatedRain[0],
    md.results.TransientSolution[-1].SmbAccumulatedRefreeze[0],
    md.results.TransientSolution[199].SmbRunoff[0],
    md.results.TransientSolution[199].SmbMelt[0],
    md.results.TransientSolution[199].SmbEC[0],
    md.results.TransientSolution[199].SmbPrecipitation[0],
    md.results.TransientSolution[199].SmbRain[0],
    md.results.TransientSolution[199].SmbRefreeze[0],
    md.results.TransientSolution[199].SmbWAdd[0]
    ]
