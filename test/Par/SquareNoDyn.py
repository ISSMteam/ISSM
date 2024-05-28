import os.path
import numpy as np
import inspect
from verbose import verbose
from transient import transient
from paterson import paterson
from arch import *

#Start defining model parameters here

#Geometry
md.geometry.thickness = 1000.0 * np.ones((md.mesh.numberofvertices))
md.geometry.base = np.zeros((md.mesh.numberofvertices))
md.geometry.surface = md.geometry.base + md.geometry.thickness

# #Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#Friction
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices))
md.friction.coefficient[np.where(md.mask.ocean_levelset < 0.)[0]] = 0.
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

#Some necessary fields to fool checkonsistency
md.initialization.vx = np.zeros((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.vz = np.zeros((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))

md.stressbalance.spcvx = np.zeros((md.mesh.numberofvertices))
md.stressbalance.spcvy = np.zeros((md.mesh.numberofvertices))
md.stressbalance.spcvz = np.zeros((md.mesh.numberofvertices))

md.stressbalance.referential = float('nan') * np.ones((md.mesh.numberofvertices, 6))
md.stressbalance.loadingforce = np.zeros((md.mesh.numberofvertices, 3))

md.smb.mass_balance = np.zeros((md.mesh.numberofvertices))

md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices))

#Numerical parameters
md.verbose = verbose(0)
md.settings.waitonlock = 30
md.groundingline.migration = 'None'

md.transient = transient.deactivateall(md.transient)

#Change name so that no test have the same name
if len(inspect.stack()) > 2:
    md.miscellaneous.name = os.path.basename(inspect.stack()[2][1]).split('.')[0]
