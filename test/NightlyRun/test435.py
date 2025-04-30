#Test Name: MISMIP3DHO
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.initialization.vx[:] = 1.
md.initialization.vy[:] = 1.
md.geometry.thickness[:] = 500. - md.mesh.x / 10000.
md.geometry.bed = -100. - md.mesh.x / 1000.
md.geometry.base = -md.geometry.thickness * md.materials.rho_ice / md.materials.rho_water
md.mask.ocean_levelset = md.geometry.thickness + md.materials.rho_water / md.materials.rho_ice * md.geometry.bed
pos = np.where(md.mask.ocean_levelset >= 0)
md.geometry.base[pos] = md.geometry.bed[pos]
md.geometry.surface = md.geometry.base + md.geometry.thickness
md = md.extrude(4, 1.)
md = setflowequation(md, 'HO', 'all')

#Boundary conditions:
md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices, ))
md.mask.ice_levelset[np.where(md.mesh.x == max(md.mesh.x))] = 0.
md.stressbalance.spcvx[:] = float('Nan')
md.stressbalance.spcvy[:] = float('Nan')
md.stressbalance.spcvz[:] = float('Nan')
posA = np.intersect1d(np.array(np.where(md.mesh.y < 1000000.1)), np.array(np.where(md.mesh.y > 999999.9)))
posB = np.intersect1d(np.array(np.where(md.mesh.y < 0.1)), np.array(np.where(md.mesh.y > -0.1)))
pos = np.unique(np.concatenate((posA, posB)))
md.stressbalance.spcvy[pos] = 0.
pos2 = np.intersect1d(np.array(np.where(md.mesh.x < 0.1)), np.array(np.where(md.mesh.x > -0.1)))
md.stressbalance.spcvx[pos2] = 0.
md.stressbalance.spcvy[pos2] = 0.

md.materials.rheology_B = 1. / ((10**-25)**(1. / 3.)) * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_law = 'None'
md.friction.coefficient[:] = np.sqrt(1e7) * np.ones((md.mesh.numberofvertices, ))
md.friction.p = 3. * np.ones((md.mesh.numberofelements, ))
md.smb.mass_balance[:] = 1.
md.basalforcings.groundedice_melting_rate[:] = 0.
md.basalforcings.floatingice_melting_rate[:] = 30.
md.transient.isthermal = 0
md.transient.isstressbalance = 1
md.transient.isgroundingline = 1
md.transient.ismasstransport = 1
md.transient.issmb = 1
md.transient.requested_outputs = ['default', 'BasalforcingsFloatingiceMeltingRate']
md.groundingline.migration = 'SubelementMigration'
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.timestepping.final_time = 30
md.timestepping.time_step = 10

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1', 'Vx1', 'Vy1', 'Vz1', 'Pressure1', 'FloatingiceMeltingrate1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2', 'Vx2', 'Vy2', 'Vz2', 'Pressure2', 'FloatingiceMeltingrate2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3', 'Vx3', 'Vy3', 'Vz3', 'Pressure3', 'FloatingiceMeltingrate3']
field_tolerances = [2e-11, 5e-12, 2e-11, 1e-11, 7e-10, 3e-08, 6e-10, 1e-13, 1e-13,
                    3e-11, 3e-11, 9e-10, 7e-11, 1e-08, 2e-07, 1e-09, 1e-10, 1e-13,
                    1e-9, 2e-08, 7e-09, 2e-7, 1e-03, 8e-04, 2e-09, 1e-10, 1e-13]
field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vz,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vz,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].BasalforcingsFloatingiceMeltingRate,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset,
                md.results.TransientSolution[2].Vx,
                md.results.TransientSolution[2].Vy,
                md.results.TransientSolution[2].Vz,
                md.results.TransientSolution[2].Pressure,
                md.results.TransientSolution[2].BasalforcingsFloatingiceMeltingRate]
