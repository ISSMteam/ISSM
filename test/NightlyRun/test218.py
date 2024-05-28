#Test Name: SquareShelfConstrainedDakotaB
import numpy as np
from model import *
from socket import gethostname
from setmask import *
from parameterize import *
from setflowequation import *
from paterson import *
from solve import *
from generic import generic
from squaremesh import *
from partitioner import *
from IssmConfig import *
from importancefactors import *

from normal_uncertain import *
from response_function import *
from dakota_method import *

from helpers import *

from verbose import *

md = squaremesh(model(), 1000000, 1000000, 5, 5)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf2.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#redo the parameter file for this special shelf.
#constant thickness, constrained (vy = 0) flow into an icefront,
#from 0 m / yr at the grounding line.

#needed later
ymin = min(md.mesh.y)
ymax = max(md.mesh.y)
xmin = min(md.mesh.x)
xmax = max(md.mesh.x)

di = md.materials.rho_ice / md.materials.rho_water

h = 1000.
md.geometry.thickness = h * np.ones((md.mesh.numberofvertices, ))
md.geometry.base = -md.materials.rho_ice / md.materials.rho_water * md.geometry.thickness
md.geometry.surface = md.geometry.base + md.geometry.thickness

#Initial velocity and pressure
md.initialization.vx = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vy = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vz = np.zeros((md.mesh.numberofvertices, ))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices, ))

#Materials
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices, ))
md.materials.rheology_B = paterson(md.initialization.temperature)
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements, ))

#Boundary conditions:
md.stressbalance.spcvx = float('Nan') * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvy = float('Nan') * np.ones((md.mesh.numberofvertices, ))
md.stressbalance.spcvz = float('Nan') * np.ones((md.mesh.numberofvertices, ))

#constrain flanks to 0 normal velocity
pos = np.where((md.mesh.x == xmin) | (md.mesh.x == xmax))
md.stressbalance.spcvx[pos] = 0.
md.stressbalance.spcvz[pos] = float('Nan')

#constrain grounding line to 0 velocity
pos = np.where(md.mesh.y == ymin)
md.stressbalance.spcvx[pos] = 0.
md.stressbalance.spcvy[pos] = 0.

#partitioning
npart = md.mesh.numberofvertices
partition = partitioner(md, 'package', 'linear', 'npart', npart)[0] - 1

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
# returns tuple "(u'6.2', )" -> unicode string '6.2', convert to float
version = float(version[0])

#variables
md.qmu.variables.rheology_B = normal_uncertain.normal_uncertain(
	'descriptor', 'scaled_MaterialsRheologyB',
	'mean', np.ones((npart, 1)),
	'stddev', .05 * np.ones((npart, 1)),
	'partition', partition
	)

#responses
md.qmu.responses.MaxVel = response_function.response_function('descriptor', 'MaxVel')

#method
md.qmu.method = dakota_method.dakota_method('nond_l')

#parameters
md.qmu.params.direct = True
md.qmu.params.interval_type = 'forward'

if version >= 6:
    md.qmu.params.analysis_driver = 'matlab'
    md.qmu.params.evaluation_scheduling = 'master'
    md.qmu.params.processors_per_evaluation = 2
else:
    md.qmu.params.analysis_driver = 'stressbalance'
    md.qmu.params.evaluation_concurrency = 1


#imperative!
md.stressbalance.reltol = 10**-10  #tighten for qmu analyses
md.qmu.isdakota = 1

#md.debug.valgrind = True
#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Stressbalance', 'overwrite', 'y')

#Fields and tolerances to track changes
md.qmu.results = md.results.dakota
md.results.dakota.importancefactors = importancefactors(md, 'scaled_MaterialsRheologyB', 'MaxVel', partition).reshape(-1, 1)
field_names = ['importancefactors']
field_tolerances = [1e-10]
field_values = [md.results.dakota.importancefactors]
