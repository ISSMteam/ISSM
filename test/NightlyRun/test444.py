#Test Name: SquareShelfTranForceNeg2dDakotaLocal

# TODO:
# - Figure out why test fails intermittently on Mac with "IndexError: list
# index out of range"
#

import numpy as np

from ContourToMesh import *
from dmeth_params_set import *
from socket import gethostname
from model import *
from parameterize import *
from partitioner import *
from regionaloutput import *
from setflowequation import *
from setmask import *
from solve import *
from triangle import *

#model not consistent:  equality thickness = surface-base violated

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.geometry.bed = md.geometry.base.copy()
pos = np.where(md.mask.ocean_levelset < 0)
md.geometry.bed[pos] = md.geometry.base[pos] - 10
md.friction.coefficient = 20. * np.ones((md.mesh.numberofvertices, ))
md.friction.p = np.ones((md.mesh.numberofelements, ))
md.friction.q = np.ones((md.mesh.numberofelements, ))
md.transient.isthermal = 0
md.transient.isgroundingline = 1
md.groundingline.migration = 'AggressiveMigration'

md.settings.output_frequency = 3
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

regionalmask = np.zeros((md.mesh.numberofvertices, ))
c_in = ContourToMesh(md.mesh.elements, md.mesh.x, md.mesh.y, '../Exp/SquareHalfRight.exp', 'node', 1)
regionalmask[np.where(c_in)] = 1
md.transient.requested_outputs = ['default', 'GroundedArea', 'FloatingArea', 'IceVolumeAboveFloatation', 'GroundedArea1', 'FloatingArea1', 'TotalFloatingBmb1', 'TotalGroundedBmb1', 'TotalSmb1',
                                  'IceMass1', 'IceVolume1', 'IceVolumeAboveFloatation1', 'IceVolumeAboveFloatation']
md.outputdefinition.definitions.append(regionaloutput('name', 'GroundedArea1', 'outputnamestring', 'GroundedArea', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition1'))
md.outputdefinition.definitions.append(regionaloutput('name', 'FloatingArea1', 'outputnamestring', 'FloatingArea', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition2'))
md.outputdefinition.definitions.append(regionaloutput('name', 'TotalFloatingBmb1', 'outputnamestring', 'TotalFloatingBmb', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition3'))
md.outputdefinition.definitions.append(regionaloutput('name', 'TotalGroundedBmb1', 'outputnamestring', 'TotalGroundedBmb', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition4'))
md.outputdefinition.definitions.append(regionaloutput('name', 'IceMass1', 'outputnamestring', 'IceMass', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition5'))
md.outputdefinition.definitions.append(regionaloutput('name', 'IceVolume1', 'outputnamestring', 'IceVolume', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition6'))
md.outputdefinition.definitions.append(regionaloutput('name', 'IceVolumeAboveFloatation1', 'outputnamestring', 'IceVolumeAboveFloatation', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition7'))
md.outputdefinition.definitions.append(regionaloutput('name', 'TotalSmb1', 'outputnamestring', 'TotalSmb', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition8'))
md.outputdefinition.definitions.append(regionaloutput('name', 'TotalSmb2', 'outputnamestring', 'TotalSmb', 'mask', regionalmask,
                                                      'definitionstring', 'Outputdefinition9'))

md.extrude(3, 1.)
md.collapse()

#Dakota options

#dakota version
version = IssmConfig('_DAKOTA_VERSION_')
version = float(version[0])

#partitioning
npart = 10
partition = partitioner(md, 'package', 'chaco', 'npart', npart, 'weighting', 'on')[0] - 1
md.qmu.isdakota = 1

#variables
md.qmu.variables.drag_coefficient = normal_uncertain.normal_uncertain(
    'descriptor', 'scaled_BasalforcingsFloatingiceMeltingRate',
    'mean', np.ones((npart, 1)),
    'stddev', .1 * np.ones((npart, 1)),
    'partition', partition
    )


#responses
md.qmu.responses.IceMass1 = response_function.response_function('descriptor','Outputdefinition5')
md.qmu.responses.IceVolume1 = response_function.response_function('descriptor','Outputdefinition6')
md.qmu.responses.IceVolumeAboveFloatation1 = response_function.response_function('descriptor','Outputdefinition7')
md.qmu.responses.IceVolumeAboveFloatation = response_function.response_function('descriptor','IceVolumeAboveFloatation')
md.qmu.responses.GroundedArea1 = response_function.response_function('descriptor','Outputdefinition1')
md.qmu.responses.FloatingArea1 = response_function.response_function('descriptor','Outputdefinition2')
md.qmu.responses.TotalFloatingBmb1 = response_function.response_function('descriptor','Outputdefinition3')
md.qmu.responses.TotalGroundedBmb1 = response_function.response_function('descriptor','Outputdefinition4')
md.qmu.responses.TotalSmb1 = response_function.response_function('descriptor','Outputdefinition8')
md.qmu.responses.TotalSmb2 = response_function.response_function('descriptor','Outputdefinition9')
md.qmu.responses.FloatingArea = response_function.response_function('descriptor','FloatingArea')

#method
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(md.qmu.method, 'seed', 1234, 'samples', 20, 'sample_type', 'random')

#parameters
md.qmu.params.direct = True
md.qmu.params.analysis_components = ''
md.qmu.params.tabular_graphics_data = True

if version >= 6:
    md.qmu.params.analysis_driver = 'matlab'
    md.qmu.params.evaluation_scheduling = 'master'
    md.qmu.params.processors_per_evaluation = 2
else:
    md.qmu.params.analysis_driver = 'stressbalance'
    md.qmu.params.evaluation_concurrency = 1

md.stressbalance.reltol = 10**-5  #tighten for qmu analyses

#solve
md.verbose = verbose('000000000')  # this line is recommended
md = solve(md, 'Transient', 'overwrite', 'y')

#Fields and tolerances to track changes
md.qmu.results = md.results.dakota

#we put all the mean and stdev data in the montecarlo field, which we will use to test for success.
md.results.dakota.montecarlo = []
for i in range(11):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].mean)

for i in range(11):
    md.results.dakota.montecarlo.append(md.results.dakota.dresp_out[i].stddev)

field_names = ['montecarlo']
field_tolerances = [1e-11]
field_values = [md.results.dakota.montecarlo]
