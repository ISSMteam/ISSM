#Test Name: SquareShelfTranIspddSicopolisMonthlyAnomaliesSSA2d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

from SMBpddSicopolis import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')

# Use of SMBpddSicopolis
md.smb = SMBpddSicopolis()
# initalize pdd fields
md.smb.initialize(md)
md.smb.s0p = md.geometry.surface.reshape(-1, 1)
md.smb.s0t = md.geometry.surface.reshape(-1, 1)

md.smb.monthlytemperatures = np.empty((md.mesh.numberofvertices, 12))
md.smb.precipitation = np.empty((md.mesh.numberofvertices, 12))
temp_ma_present = -10. * np.ones((md.mesh.numberofvertices, )) - md.smb.rlaps * md.geometry.surface / 1000.
temp_mj_present = 10. * np.ones((md.mesh.numberofvertices, )) - md.smb.rlaps * md.geometry.surface / 1000.
precipitation = 5. * np.ones((md.mesh.numberofvertices, ))
for imonth in range(12):
    md.smb.monthlytemperatures[0:md.mesh.numberofvertices, imonth] = md.materials.meltingpoint + temp_ma_present + (temp_mj_present - temp_ma_present) * np.sin((imonth + 1. - 4.) * np.pi / 6.0)
    md.smb.precipitation[0:md.mesh.numberofvertices, imonth] = precipitation

# time steps and resolution
md.timestepping.time_step = 1
md.settings.output_frequency = 1
md.timestepping.final_time = 2

md.transient.issmb = 1
md.transient.ismasstransport = 1
md.transient.isstressbalance = 0
md.transient.isthermal = 0

md.transient.requested_outputs = ['default', 'TemperaturePDD']
md.cluster = generic('name', gethostname(), 'np', 1)  # 3 for the cluster

nv = md.mesh.numberofvertices


def grab(md):
    """Copy the fields we track so that a later solve cannot overwrite them"""
    steps = md.results.TransientSolution
    return [(np.copy(steps[i].TemperaturePDD), np.copy(steps[i].SmbMassBalance)) for i in range(2)]

# Run 1: legacy anomalies (one value per vertex, same for all 12 months)
md.smb.temperature_anomaly = 1.5 * np.ones((nv, 1))
md.smb.precipitation_anomaly = 1.1 * np.ones((nv, 1))
md = solve(md, 'Transient')
res1 = grab(md)

# Run 2: monthly anomalies (list of 12 time series, time row in years last) holding the same constants.
# Must reproduce Run 1 exactly
time_row = np.array([[0., 2.]])
md.smb.temperature_anomaly = [np.vstack((1.5 * np.ones((nv, 2)), time_row)) for imonth in range(12)]
md.smb.precipitation_anomaly = [np.vstack((1.1 * np.ones((nv, 2)), time_row)) for imonth in range(12)]
md = solve(md, 'Transient')
res2 = grab(md)

# Run 3: anomalies that vary with calendar month and with time
tanom = []
panom = []
for imonth in range(1, 13):
    dT = 2. * np.sin(2. * np.pi * (imonth - 1) / 12.)
    tanom.append(np.vstack((np.hstack((dT * np.ones((nv, 1)), (dT + 1.) * np.ones((nv, 1)))), time_row)))
    panom.append(np.vstack((np.hstack(((1. + 0.02 * imonth) * np.ones((nv, 1)), (1.2 + 0.02 * imonth) * np.ones((nv, 1)))), time_row)))
md.smb.temperature_anomaly = tanom
md.smb.precipitation_anomaly = panom
md = solve(md, 'Transient')
res3 = grab(md)

# Consistency check: 12 identical constant monthly series must reproduce the legacy results
for i in range(2):
    for k in range(2):
        if np.max(np.abs(res1[i][k] - res2[i][k])) > 1e-12 * np.max(np.abs(res1[i][k])):
            raise Exception('monthly anomalies do not reproduce the legacy (nv x 1) anomaly results')

#Fields and tolerances to track changes
field_names = [
    'LegacyTemperaturePDD1', 'LegacySmbMassBalance1', 'LegacyTemperaturePDD2', 'LegacySmbMassBalance2',
    'MonthlyConstTemperaturePDD1', 'MonthlyConstSmbMassBalance1', 'MonthlyConstTemperaturePDD2', 'MonthlyConstSmbMassBalance2',
    'MonthlyVaryingTemperaturePDD1', 'MonthlyVaryingSmbMassBalance1', 'MonthlyVaryingTemperaturePDD2', 'MonthlyVaryingSmbMassBalance2']
field_tolerances = [1e-13] * 12
field_values = [
    res1[0][0], res1[0][1], res1[1][0], res1[1][1],
    res2[0][0], res2[0][1], res2[1][0], res2[1][1],
    res3[0][0], res3[0][1], res3[1][0], res3[1][1]]
