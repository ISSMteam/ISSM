import numpy as np
from SetIceSheetBC import SetIceSheetBC

#Ok, start defining model parameters here
md.verbose = 2

print("      creating thickness")
md.geometry.surface = -md.mesh.x * np.tan(3. * np.pi / 180.)
#md.geometry.base = md.geometry.surface-1000.
md.geometry.base = md.geometry.surface - 1000. + 100. * np.exp(-((md.mesh.x - np.max(md.mesh.x) / 2.)**2 + (md.mesh.y - np.max(md.mesh.y) / 2.)**2) / (10000.**2))
md.geometry.thickness = md.geometry.surface - md.geometry.base

print("      creating drag")
md.friction.coefficient = np.sqrt(md.constants.yts / (2.140373 * 10**-7 * 1000.)) * np.ones((md.mesh.numberofvertices))
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.zeros((md.mesh.numberofelements))

print("      creating flow law parameter")
md.materials.rheology_B = 1.4734 * 10**14 * np.ones((md.mesh.numberofvertices))
md.materials.rheology_n = 1. * np.ones((md.mesh.numberofelements))
md.materials.rheology_law = 'None'

print("      boundary conditions for stressbalance model")
#Create node on boundary first (because we cannot use mesh)
md = SetIceSheetBC(md)
md.stressbalance.spcvx = 100. * np.ones((md.mesh.numberofvertices))
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.vel = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))
md.initialization.temperature = 255. * np.ones((md.mesh.numberofvertices))
pos = np.nonzero(np.logical_or(np.logical_or(md.mesh.x == np.min(md.mesh.x), md.mesh.x == np.max(md.mesh.x)), np.logical_or(md.mesh.y == np.min(md.mesh.y), md.mesh.y == np.max(md.mesh.y))))
md.balancethickness.spcthickness = float('NaN') * np.ones((md.mesh.numberofvertices))
md.balancethickness.spcthickness[pos] = md.geometry.thickness[pos]
md.masstransport.spcthickness = float('NaN') * np.ones((md.mesh.numberofvertices))
md.masstransport.spcthickness[pos] = md.geometry.thickness[pos]
md.thermal.spctemperature = 255. * np.ones((md.mesh.numberofvertices))
md.basalforcings.geothermalflux = 0.4 * np.ones((md.mesh.numberofvertices))

#Parallel options
md.mesh.average_vertex_connectivity = 200

#Transient options
md.timestepping.time_step = 1.
md.timestepping.final_time = 10.
md.masstransport.stabilization = 1
md.thermal.stabilization = 1
md.thermal.penalty_threshold = 10**5
md.transient.isthermal = 0
