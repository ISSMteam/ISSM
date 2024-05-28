#Test Name: SquareShelfTranIspddIsdeltaO18pdNoInterpSSA2d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from cuffey import *
from frictionjosh import *

from SMBpddSicopolis import *

md = triangle(model(), '../Exp/Square.exp', 180000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')

# Use of ispdd and isdelta18o methods
md.smb = SMBd18opdd()
md.smb.isd18opd=1

# Add temperature, precipitation and delta18o needed to measure the surface mass balance
# creating delta18o
delta18o = np.loadtxt('../Data/delta18o.data')
md.smb.delta18o = delta18o
md.smb.istemperaturescaled = 0
md.smb.isprecipscaled = 0

# creating Present day temperatures
# Same temperature over the all region:
curve=np.sin(np.linspace(0, np.pi, 12))
tmonth = np.ones(12) * (238.15 + 20.) + 20.*curve
md.smb.temperatures_presentday = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.temperatures_presentday[0:md.mesh.numberofvertices, imonth] = tmonth[imonth]

md.smb.initialize(md)

# creating precipitation
md.smb.precipitations_presentday = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.precipitations_presentday[0:md.mesh.numberofvertices, imonth] = -0.4 * 10**(-6) * md.mesh.y + 0.5

# 3 total years of input
md.smb.temperatures_reconstructed = np.empty((md.mesh.numberofvertices+1, 12*3))
md.smb.precipitations_reconstructed = np.empty((md.mesh.numberofvertices+1, 12*3))

md.smb.temperatures_reconstructed[0:md.mesh.numberofvertices,0:12] = md.smb.temperatures_presentday
md.smb.temperatures_reconstructed[0:md.mesh.numberofvertices,12:24] = md.smb.temperatures_presentday+1.2
md.smb.temperatures_reconstructed[0:md.mesh.numberofvertices,24:36] = md.smb.temperatures_presentday-0.8

md.smb.precipitations_reconstructed[0:md.mesh.numberofvertices,0:12] = md.smb.precipitations_presentday
md.smb.precipitations_reconstructed[0:md.mesh.numberofvertices,12:24] = md.smb.precipitations_presentday+0.1
md.smb.precipitations_reconstructed[0:md.mesh.numberofvertices,24:36] = md.smb.precipitations_presentday-0.1

tim1 = np.linspace(1,12,12)/12

md.smb.temperatures_reconstructed[md.mesh.numberofvertices,0:12] = tim1
md.smb.temperatures_reconstructed[md.mesh.numberofvertices,12:24] = tim1+3
md.smb.temperatures_reconstructed[md.mesh.numberofvertices,24:36] = tim1+5

md.smb.precipitations_reconstructed[md.mesh.numberofvertices,0:12] = tim1
md.smb.precipitations_reconstructed[md.mesh.numberofvertices,12:24] = tim1+3
md.smb.precipitations_reconstructed[md.mesh.numberofvertices,24:36] = tim1+5

# creating initialization and spc temperatures initialization and
# spc
md.thermal.spctemperature=np.mean(md.smb.temperatures_presentday[0:md.mesh.numberofvertices,0:12],axis=1)-10.
md.initialization.temperature=md.thermal.spctemperature

md.smb.s0p = np.maximum(md.geometry.surface.reshape(-1, 1),np.zeros((md.mesh.numberofvertices,1)))
md.smb.s0t = np.maximum(md.geometry.surface.reshape(-1, 1),np.zeros((md.mesh.numberofvertices,1)))
md.smb.issetpddfac = 1
md.smb.pddfac_snow = 8
md.smb.pddfac_ice = 10

md.extrude(5,1.2)
md=setflowequation(md,'HO','all')
md.settings.results_on_nodes=['Temperature','Waterfraction','Enthalpy']

md.thermal.isenthalpy=1
md.thermal.isdynamicbasalspc=1
md.thermal.fe = 'P1xP2'

md.initialization.waterfraction = np.zeros((md.mesh.numberofvertices,1))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices,1))
md.initialization.enthalpy = np.zeros((md.mesh.numberofvertices,1))
md.thermal.isdrainicecolumn = 0

md = solve(md, 'thermal')

md.initialization.temperature = md.results.ThermalSolution.Temperature
md.initialization.enthalpy = md.results.ThermalSolution.Enthalpy
md.materials.rheology_B = cuffey(md.initialization.temperature).reshape(np.shape(md.initialization.temperature))

# Friction
TEMP = np.zeros((md.mesh.numberofvertices,))
TEMP[md.mesh.elements - 1] = md.initialization.temperature[:,0:6].reshape(md.mesh.numberofelements,6)

temperature = TEMP
pressure = md.initialization.pressure
Tm = md.materials.meltingpoint-md.materials.beta*pressure

coefficient=md.friction.coefficient
md.friction=frictionjosh()
md.friction.coefficient = coefficient
md.friction.pressure_adjusted_temperature = temperature - Tm
md.friction.gamma = 5

# time steps and resolution
md.timestepping.time_step = 0.5
md.settings.output_frequency = 1
md.timestepping.final_time = 2
md.timestepping.interp_forcing = 0

md.transient.requested_outputs = ['default', 'IceVolumeAboveFloatation','IceVolume','TemperaturePDD']
md=setflowequation(md,'SSA','all')
md.cluster = generic('name', gethostname(), 'np', 1)  # 3 for the cluster
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'Temperature1', 'Enthalpy1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'Temperature2', 'Enthalpy2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'Temperature3', 'Enthalpy3', 'SmbMassBalance3',
               'Vx4', 'Vy4', 'Vz4', 'Vel4', 'Pressure4', 'Bed4', 'Surface4', 'Thickness4', 'Temperature4', 'Enthalpy4', 'SmbMassBalance4']
field_tolerances=[1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
   1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
   1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
   1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].Enthalpy,
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
                md.results.TransientSolution[1].Enthalpy,
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
                md.results.TransientSolution[2].Enthalpy,
                md.results.TransientSolution[2].SmbMassBalance,
                md.results.TransientSolution[3].Vx,
                md.results.TransientSolution[3].Vy,
                md.results.TransientSolution[3].Vz,
                md.results.TransientSolution[3].Vel,
                md.results.TransientSolution[3].Pressure,
                md.results.TransientSolution[3].Base,
                md.results.TransientSolution[3].Surface,
                md.results.TransientSolution[3].Thickness,
                md.results.TransientSolution[3].Temperature,
                md.results.TransientSolution[3].Enthalpy,
                md.results.TransientSolution[3].SmbMassBalance]
