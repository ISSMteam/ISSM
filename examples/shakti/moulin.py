import numpy as np
from SetIceSheetBC import SetIceSheetBC
from frictionshakti import frictionshakti
from verbose import verbose

#Start defining model parameters here
# Set up bed topography and ice geometry for a tilted 500m thick slab
md.geometry.base = .02 * md.mesh.x
md.geometry.bed = md.geometry.base
md.geometry.surface = .02 * md.mesh.x + 500
md.geometry.thickness = md.geometry.surface - md.geometry.bed

# Define ice sliding velocity (m / yr)
md.initialization.vx = 1.0e-6 * md.constants.yts * np.ones(md.mesh.numberofvertices)
md.initialization.vy = np.zeros(md.mesh.numberofvertices)
md.initialization.pressure = np.zeros(md.mesh.numberofvertices)

# Materials
# Ice flow law parameter (note that the standard parameter A = B^(- 3))
md.materials.rheology_B = 5e-25**(-1. / 3.) * np.ones(md.mesh.numberofvertices)
md.initialization.temperature = 273. * np.ones(md.mesh.numberofvertices)
md.materials.rheology_n = 3. * np.ones(md.mesh.numberofelements)

#Calving
md.calving.calvingrate = np.zeros(md.mesh.numberofvertices)
#md.calving.spclevelset = NaN(md.mesh.numberofvertices, 1)

# Friction law and coefficient
md.friction = frictionshakti(md)
md.friction.coefficient = 20. * np.ones(md.mesh.numberofvertices)

#Numerical parameters
#md.stressbalance.viscosity_overshoot = 0.0
md.masstransport.stabilization = 1.
md.thermal.stabilization = 1.
md.verbose = verbose(0)
md.settings.waitonlock = 30
md.stressbalance.restol = 0.05
md.steadystate.reltol = 0.05
md.stressbalance.reltol = 0.05
md.stressbalance.abstol = np.nan
md.timestepping.time_step = 1.
md.timestepping.final_time = 3.

#GIA:
md.gia.lithosphere_thickness = 100. * np.ones(md.mesh.numberofvertices)  # in km
md.gia.mantle_viscosity = 1.0e21 * np.ones(md.mesh.numberofvertices)  # in Pa.s
md.materials.lithosphere_shear_modulus = 6.7e10   # in Pa
md.materials.lithosphere_density = 3.32   # in g / cm^ - 3
md.materials.mantle_shear_modulus = 1.45e11  # in Pa
md.materials.mantle_density = 3.34      # in g / cm^ - 3

#Boundary conditions:
md = SetIceSheetBC(md)

#Change name so that no test have the same name
md.private.runtimename = True
