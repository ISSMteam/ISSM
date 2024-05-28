#Test Name: SquareShelfTranIspddIsdeltaSSA3d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

from generic import generic

md = triangle(model(), '../Exp/Square.exp', 600000)  #180000
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')

# Use of ispdd and isdelta18o methods
md.smb = SMBpdd()
md.smb.isdelta18o = 0
md.smb.ismungsm = 1

# time steps and resolution
md.timestepping.time_step = 20.
md.settings.output_frequency = 1
md.timestepping.final_time = 60.

# creating Present day and lgm temperatures
# Same temperature over the all region:
curve=np.sin(np.linspace(0, np.pi, 12))
tmonth = np.ones(12) * (238.15 + 20.) + 20.*curve
md.smb.temperatures_presentday = np.zeros((md.mesh.numberofvertices, 12))
md.smb.temperatures_lgm = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.temperatures_presentday[0:md.mesh.numberofvertices, imonth] = tmonth[imonth]
    md.smb.temperatures_lgm[0:md.mesh.numberofvertices, imonth] = tmonth[imonth] - 20.

# creating initialization and spc temperatures initialization and spc
md.thermal.spctemperature = np.mean(md.smb.temperatures_lgm[0:md.mesh.numberofvertices, :], axis=1)  # - 10 * ones(md.mesh.numberofvertices, 1)
md.thermal.spctemperature = np.tile(md.thermal.spctemperature, (int(md.timestepping.final_time / md.timestepping.time_step), 1)).T
itemp = np.arange(0, md.timestepping.final_time, md.timestepping.time_step)
md.thermal.spctemperature = np.vstack((md.thermal.spctemperature, itemp))

md.initialization.temperature = md.smb.temperatures_lgm[0:md.mesh.numberofvertices, 0]  # * ones(md.mesh.numberofvertices, 1)
md.smb.initialize(md)

# creating precipitation
md.smb.precipitations_presentday = np.zeros((md.mesh.numberofvertices, 12))
md.smb.precipitations_lgm = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.precipitations_presentday[0:md.mesh.numberofvertices, imonth] = -0.4 * 10**(-6) * md.mesh.y + 0.5
    md.smb.precipitations_lgm[0:md.mesh.numberofvertices, imonth] = -0.4 * 10**(-6) * md.mesh.y + 0.5

fsize = int(md.timestepping.final_time / md.timestepping.time_step) + 2
md.smb.Pfac = np.zeros((2, fsize))
md.smb.Tdiff = np.zeros((2, fsize))
md.smb.sealev = np.zeros((2, fsize))
for iint in range(0, fsize):
    # Interpolation factors
    md.smb.Pfac[0, iint] = 0.15 * (iint + 1)
    md.smb.Tdiff[0, iint] = 0.15 * (iint + 1)
    md.smb.sealev[0, iint] = 0.15 * (iint + 1)
    # Year of each data point
    md.smb.Pfac[1, iint] = (float(iint)) * 20
    md.smb.Tdiff[1, iint] = (float(iint)) * 20
    md.smb.sealev[1, iint] = (float(iint)) * 20

md.smb.issetpddfac=1.
md.smb.pddfac_snow=2.
md.smb.pddfac_ice=2.

md.extrude(3, 1.)
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'Temperature1', 'BasalforcingsGroundediceMeltingRate1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'Temperature2', 'BasalforcingsGroundediceMeltingRate2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'Temperature3', 'BasalforcingsGroundediceMeltingRate3', 'SmbMassBalance3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-8, 1e-8, 1e-8, 1e-13, 1e-8, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-8, 1e-8, 1e-8, 7e-13, 2e-7, 1e-13,
                    1e-13, 1e-13, 1e-08, 1e-13, 1e-13, 1e-8, 1e-8, 1e-8, 7e-13, 5e-7, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].BasalforcingsGroundediceMeltingRate,
                md.results.TransientSolution[2].SmbMassBalance]
