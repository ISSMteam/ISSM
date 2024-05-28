#Test Name: SquareSheetConstrainedSmbGradients3d
import copy
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md.extrude(3, 1)
md = setflowequation(md, 'HO', 'all')
md.smb = SMBgradients()
md.smb.b_pos = (-100. + 0.00005 * md.mesh.x - 0.0001 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice
md.smb.b_neg = (250. + 0.000051 * md.mesh.x - 0.00011 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice
md.smb.href = copy.deepcopy(md.geometry.surface)
md.smb.smbref = (1000. - 0.001 * md.mesh.x - 0.005 * md.mesh.y) / 1000. * md.materials.rho_freshwater / md.materials.rho_ice
md.transient.requested_outputs = ['default', 'TotalSmb']
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1',
               'Bed1', 'Surface1', 'Thickness1',
               'Temperature1', 'SMB1', 'TotalSmb1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2',
               'Bed2', 'Surface2', 'Thickness2',
               'Temperature2', 'SMB2', 'TotalSmb2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3',
               'Bed3', 'Surface3', 'Thickness3',
               'Temperature3', 'SMB3', 'TotalSmb3']
field_tolerances = [1e-09, 1e-09, 1e-09, 1e-09,
                    1e-10, 1e-10, 1e-10,
                    1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 2e-10, 1e-09,
                    1e-10, 1e-10, 1e-10,
                    1e-10, 1e-10, 1e-10,
                    1e-09, 5e-09, 1e-09, 1e-09,
                    1e-10, 1e-10, 1e-10,
                    1e-10, 1e-10, 1e-10]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].Temperature,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[0].TotalSmb,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].Temperature,
                md.results.TransientSolution[1].SmbMassBalance,
                md.results.TransientSolution[1].TotalSmb,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Vel,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].Temperature,
                md.results.TransientSolution[2].SmbMassBalance,
                md.results.TransientSolution[2].TotalSmb]
