#Test Name: SquareSheetConstrainedSmbGradientsEla3d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from SMBgradientsela import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')

#Change geometry
md.geometry.thickness = md.geometry.surface * 30.
md.geometry.surface = md.geometry.base + md.geometry.thickness

md = md.extrude(3, 1.)
md = setflowequation(md, 'HO', 'all')
md.smb = SMBgradientsela()
md.smb.ela = 1500. * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_pos = 0.002 * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_neg = 0.005 * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_max = 4. * (md.materials.rho_freshwater / md.materials.rho_ice) * np.ones((md.mesh.numberofvertices + 1, ))
md.smb.b_min = -4. * (md.materials.rho_freshwater / md.materials.rho_ice) * np.ones((md.mesh.numberofvertices + 1, ))


#Transient options
md.transient.requested_outputs = ['default', 'TotalSmb']
md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vz1', 'Vel1', 'Bed1', 'Surface1', 'Thickness1', 'Temperature1', 'SMB1', 'TotalSmb1',
               'Vx2', 'Vy2', 'Vz2', 'Vel2', 'Bed2', 'Surface2', 'Thickness2', 'Temperature2', 'SMB2', 'TotalSmb2',
               'Vx3', 'Vy3', 'Vz3', 'Vel3', 'Bed3', 'Surface3', 'Thickness3', 'Temperature3', 'SMB3', 'TotalSmb3']
field_tolerances = [1e-09, 1e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 1e-09, 1e-10, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10,
                    1e-09, 5e-09, 1e-09, 1e-09, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]
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
