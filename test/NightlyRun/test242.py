#Test Name: SquareShelfTranForceNoInterp3d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


md = triangle(model(), '../Exp/Square.exp', 350000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')
md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

md.timestepping.time_step = 1.
md.settings.output_frequency = 1
md.timestepping.final_time = 4.
md.timestepping.interp_forcing = False

#Set up transient
smb = np.ones((md.mesh.numberofvertices)) * 3.6
smb = np.vstack((smb, smb * -1.)).T

md.smb.mass_balance = np.vstack((smb, [1.5, 3.]))
md.transient.isthermal = False

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Pressure1', 'Bed1', 'Surface1', 'Thickness1', 'SmbMassBalance1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Pressure2', 'Bed2', 'Surface2', 'Thickness2', 'SmbMassBalance2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Pressure3', 'Bed3', 'Surface3', 'Thickness3', 'SmbMassBalance3',
               'Vx4', 'Vy4', 'Vz4', 'Vel4', 'Pressure4', 'Bed4', 'Surface4', 'Thickness4', 'SmbMassbalance4']
field_tolerances = [1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
                    1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
                    1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13,
                    1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].SmbMassBalance,
                md.results.TransientSolution[3].Vx,
                md.results.TransientSolution[3].Vy,
                md.results.TransientSolution[3].Vz,
                md.results.TransientSolution[3].Vel,
                md.results.TransientSolution[3].Pressure,
                md.results.TransientSolution[3].Base,
                md.results.TransientSolution[3].Surface,
                md.results.TransientSolution[3].Thickness,
                md.results.TransientSolution[3].SmbMassBalance]
