#Test Name: SquareSheetShelfDakotaScaledResponse
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from partitioner import *

md = triangle(model(), '../Exp/Square.exp', 200000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#partitioning
npart = 10
partition = partitioner(md, 'package', 'chaco', 'npart', npart)[0] - 1
md.qmu.isdakota = 1

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#variables
md.qmu.variables.rho_ice = normal_uncertain.normal_uncertain(
    'descriptor', 'MaterialsRhoIce',
    'mean', 1,
    'stddev', .01
    )

#responses
md.qmu.responses.MaxVel = response_function.response_function(
    'descriptor', 'scaled_Thickness',
    'partition', partition
    )

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
md.stressbalance.reltol = 10**-5  #tighten for qmu analysese

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Stressbalance', 'overwrite', 'y')
md.qmu.results = md.results.dakota

#test on thickness
h = np.zeros(npart)
for i in range(npart):
    h[i] = md.qmu.results.dresp_out[i].mean

#project onto grid
thickness = h[partition]

#Fields and tolerances to track changes
field_names = ['Thickness']
field_tolerances = [1e-10]
field_values = [thickness]
