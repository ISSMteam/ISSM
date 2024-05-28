#Test Name: SquareSheetConstrainedSmbMeltComponentsClim2d
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
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

md.timestepping.time_step = 0.5
md.settings.output_frequency = 1
md.timestepping.final_time = 6.

#Set up transient
smb = np.ones((md.mesh.numberofvertices)) * 3.6
smb = np.vstack((smb, smb * -1.)).T

md.smb = SMBmeltcomponents()
md.smb.accumulation = np.vstack((smb, [1.5, 3.]))
md.smb.evaporation = np.vstack((smb / 2, [1.5, 3.]))
md.smb.melt = np.vstack((smb / 2, [1.5, 3.]))
md.smb.refreeze = np.vstack((smb, [1.5, 3.]))
md.transient.isthermal = False

md.timestepping.cycle_forcing = 1

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'SmbMassBalance3',
               'Vx4', 'Vy4', 'Vel4', 'Pressure4', 'Bed4', 'Surface4', 'Thickness4', 'SmbMassBalance4']
field_tolerances = [1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[3].Vx,
                md.results.TransientSolution[3].Vy,
                md.results.TransientSolution[3].Vel,
                md.results.TransientSolution[3].Pressure,
                md.results.TransientSolution[3].Base,
                md.results.TransientSolution[3].Surface,
                md.results.TransientSolution[3].Thickness,
                md.results.TransientSolution[3].SmbMassBalance,
                md.results.TransientSolution[7].Vx,
                md.results.TransientSolution[7].Vy,
                md.results.TransientSolution[7].Vel,
                md.results.TransientSolution[7].Pressure,
                md.results.TransientSolution[7].Base,
                md.results.TransientSolution[7].Surface,
                md.results.TransientSolution[7].Thickness,
                md.results.TransientSolution[7].SmbMassBalance,
                md.results.TransientSolution[11].Vx,
                md.results.TransientSolution[11].Vy,
                md.results.TransientSolution[11].Vel,
                md.results.TransientSolution[11].Pressure,
                md.results.TransientSolution[11].Base,
                md.results.TransientSolution[11].Surface,
                md.results.TransientSolution[11].Thickness,
                md.results.TransientSolution[11].SmbMassBalance]
