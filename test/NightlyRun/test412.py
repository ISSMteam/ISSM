#Test Name: SquareSheetShelfDiadSSA3dDakota
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from partitioner import *
from importancefactors import *
from normal_uncertain import *
from response_function import *

md = triangle(model(), '../Exp/Square.exp', 300000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#partitioning
npart = md.mesh.numberofvertices
partition = partitioner(md, 'package', 'linear', 'npart', npart)[0] - 1
md.qmu.isdakota = 1

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#variables
md.qmu.variables.rho_ice = normal_uncertain.normal_uncertain(
    'descriptor', 'MaterialsRhoIce',
    'mean', md.materials.rho_ice,
    'stddev', .01
    )
md.qmu.variables.drag_coefficient = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_FrictionCoefficient',
    'mean', np.ones((npart, 1)),
    'stddev', .01 * np.ones((npart, 1)),
    'partition', partition
    )

#responses
md.qmu.responses.MaxVel = response_function.response_function('descriptor','MaxVel')

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
md.stressbalance.reltol = 10**-5  #tighten for qmu analyses

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Stressbalance', 'overwrite', 'y')

#Fields and tolerances to track changes
md.qmu.results = md.results.dakota
md.results.dakota.importancefactors = importancefactors(md, 'scaled_FrictionCoefficient', 'MaxVel', partition).T
field_names = ['importancefactors']
field_tolerances = [1e-10]
field_values = [md.results.dakota.importancefactors]
