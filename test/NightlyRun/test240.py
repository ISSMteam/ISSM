#Test Name: SquareShelfTranIspddIsdeltaO18pdNoInterpSSA2d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')

# Use of ispdd and isdelta18o methods
md.smb = SMBd18opdd()
md.smb.isd18opd = 1

# Add temperature, precipitation and delta18o needed to measure the surface mass balance
# creating delta18o
delta18o = np.loadtxt('../Data/delta18o.data')
md.smb.delta18o = delta18o

# creating Present day temperatures
# Same temperature over the all region:
tmonth = np.ones(12) * (238.15 + 20.)
md.smb.temperatures_presentday = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.temperatures_presentday[0:md.mesh.numberofvertices, imonth] = tmonth[imonth]

# creating initialization and spc temperatures initialization and spc
md.thermal.spctemperature = np.mean(md.smb.temperatures_presentday[0:md.mesh.numberofvertices, :], axis=1).reshape(-1, 1)
md.thermal.spctemperature = md.thermal.spctemperature - 10
md.initialization.temperature = md.thermal.spctemperature
md.smb.initialize(md)

# creating precipitation
md.smb.precipitations_presentday = np.zeros((md.mesh.numberofvertices, 12))
for imonth in range(0, 12):
    md.smb.precipitations_presentday[0:md.mesh.numberofvertices, imonth] = -0.4 * 10**(-6) * md.mesh.y + 0.5

# time steps and resolution
md.timestepping.time_step = 0.5
md.settings.output_frequency = 1
md.timestepping.final_time = 2
md.timestepping.interp_forcing = 0

md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'SmbMassBalance3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].SmbMassBalance]
