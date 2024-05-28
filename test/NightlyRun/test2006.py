#Test Name: EarthSlc Dakota Sampling glaciers
import numpy as np
import pickle
from socket import gethostname
from dmeth_params_set import *
from gmtmask import *
from lovenumbers import *
from materials import *
from MatlabFuncs import *
from model import *
from nodalvalue import *
from normal_uncertain import *
from response_function import *
from solve import *


# Mesh earth
md = model()
md.cluster = generic('name', gethostname(), 'np', 5)
md.mesh = gmshplanet('radius', 6.371012 * 1e3, 'resolution', 700.)  #700 km resolution mesh

# Load precomputed mesh
with open('../Data/SlcTestMesh.pkl', 'rb') as slc_test_mesh_file:
    md.mesh = pickle.load(slc_test_mesh_file)

# Geometry for the bed, arbitrary thickness of 100
md.geometry.bed = np.zeros((md.mesh.numberofvertices, ))
md.geometry.base = md.geometry.bed
md.geometry.thickness = 100 * np.ones((md.mesh.numberofvertices, ))
md.geometry.surface = md.geometry.bed + md.geometry.thickness

# Solidearth loading #{{{
md.masstransport.spcthickness = np.append(md.geometry.thickness, 0)
md.smb.mass_balance = np.zeros((md.mesh.numberofvertices, ))
# Antarctica
xe = md.mesh.x[md.mesh.elements - 1].sum(axis=1) / 3
ye = md.mesh.y[md.mesh.elements - 1].sum(axis=1) / 3
ze = md.mesh.z[md.mesh.elements - 1].sum(axis=1) / 3
re = pow((pow(xe, 2) + pow(ye, 2) + pow(ze, 2)), 0.5)

late = asind(ze / re)
longe = atan2d(ye, xe)
pos = np.where(late < -80)[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 100
posant = pos
# Greenland
pos = np.where(np.logical_and.reduce((late > 60, late < 90, longe > -75, longe < -15)))[0]
md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] = md.masstransport.spcthickness[md.mesh.elements[pos, :] - 1] - 100
posgre = pos

# Elastic loading from love numbers
md.solidearth.lovenumbers = lovenumbers('maxdeg', 100)
#}}}

# Mask #{{{
mask = gmtmask(md.mesh.lat, md.mesh.long)
oceanmask = -1 * np.ones((md.mesh.numberofvertices, 1))
pos = np.where(mask == 0)[0]
oceanmask[pos] = 1

icemask = np.ones((md.mesh.numberofvertices, 1))
# NOTE: Need to be careful here: when addressing with multidimensional array in
# MATLAB, only first column of values are used as indices
#
icemask[md.mesh.elements[posant][:,0] - 1] = -1
icemask[md.mesh.elements[posgre][:,0] - 1] = -1

md.mask.ice_levelset = icemask
md.mask.ocean_levelset = oceanmask

# Time stepping
md.timestepping.start_time = 0
md.timestepping.time_step = 1
md.timestepping.final_time = 10

# Masstransport
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vx = np.zeros((md.mesh.numberofvertices, ))
md.initialization.vy = np.zeros((md.mesh.numberofvertices, ))
md.initialization.sealevel = np.zeros((md.mesh.numberofvertices, ))
md.initialization.str = 0

# Materials
md.materials = materials('hydro')

# Miscellaneous
md.miscellaneous.name = 'test2006'

# Solution parameters
md.cluster.np = 3
md.solidearth.settings.reltol = np.nan
md.solidearth.settings.abstol = 1e-3
md.solidearth.settings.sealevelloading = 1
md.solidearth.settings.isgrd = 1
md.solidearth.settings.ocean_area_scaling = 0
md.solidearth.settings.grdmodel = 1

md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 1
md.solidearth.settings.viscous = 0

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1
md.solidearth.requested_outputs = ['Sealevel']

dh = md.masstransport.spcthickness
deltathickness = np.zeros((md.mesh.numberofvertices + 1, 10 + 1)) # NOTE: Need to add another row as it is added in loop in MATLAB, which causes a RuntimeError in Python
for i in range(10 + 1):
    deltathickness[0:-1, i] = md.geometry.thickness + dh[0:-1] * i
deltathickness[-1, :] = np.arange(0, 10 + 1, 1)
md.masstransport.spcthickness = deltathickness

# Hack
md.geometry.surface = np.zeros((md.mesh.numberofvertices, 1))
md.geometry.thickness = np.ones((md.mesh.numberofvertices, 1))
md.geometry.base = -np.ones((md.mesh.numberofvertices, 1))
md.geometry.bed = md.geometry.base

# Uncertainty quantification
# Ice sheets #{{{
npart = 1
nt = 1
partition = -1 * np.ones((md.mesh.numberofelements, 1))
pos = np.where(late < -80)[0]
partition[pos] = 0
pos = np.where(np.logical_and.reduce((late > 70, late < 80, longe > -60, longe < -30)))[0]
partition[pos] = 0

# Variables
qmuvar = OrderedStruct()
qmuvar.surfaceload = normal_uncertain.normal_uncertain(
    'descriptor',   'scaled_SurfaceloadIceThicknessChange',
    'mean',         1 * np.ones((npart, nt)),
    'stddev',       1 * np.ones((npart, nt)), # 10% standard deviation
    'partition',    partition,
    'transient',    'on',
    'nsteps',       nt
)
#}}}

# Correlation
md.qmu.correlation_matrix = []

# Variables final declaration
md.qmu.variables = OrderedStruct()
md.qmu.variables.surfaceload = qmuvar.surfaceload

locations = np.array([1, 5, 10, 15, 20])

# Responses #{{{
md.qmu.responses.sealevel1 = response_function.response_function('descriptor', 'Outputdefinition1')
md.qmu.responses.sealevel2 = response_function.response_function('descriptor', 'Outputdefinition2')
md.qmu.responses.sealevel3 = response_function.response_function('descriptor', 'Outputdefinition3')
md.qmu.responses.sealevel4 = response_function.response_function('descriptor', 'Outputdefinition4')
md.qmu.responses.sealevel5 = response_function.response_function('descriptor', 'Outputdefinition5')

# Output definitions
for i in range(len(locations)):
    if i == 0:
        md.outputdefinition.definitions = [nodalvalue('name', 'SNode',
                                                      'definitionstring', 'Outputdefinition1',
                                                      'model_string', 'Sealevel',
                                                      'node', locations[i])]
    else:
        md.outputdefinition.definitions.append(
            nodalvalue('name', 'SNode',
                       'definitionstring', 'Outputdefinition' + str(i + 1),
                       'model_string', 'Sealevel',
                       'node', locations[i]))
#}}}

# Algorithm #{{{
md.qmu.method = dakota_method.dakota_method('nond_samp')
md.qmu.method = dmeth_params_set(md.qmu.method,
                                 'seed', 1234,
                                 'samples', 10,
                                 'sample_type', 'random')
md.qmu.output = 1
#}}}
# Parameters #{{{
md.qmu.params.direct = True
md.qmu.params.interval_type = 'forward'
md.qmu.params.analysis_driver = 'matlab'
md.qmu.params.evaluation_scheduling = 'master'
md.qmu.params.processors_per_evaluation = 2
md.qmu.params.tabular_graphics_data = True
md.qmu.isdakota = 1
md.verbose = verbose(0)
md.verbose.qmu = 1
#}}}
# QMU statistics #{{{

# TODO: Abstract reshaping of arrays away to src/m/classes/qmustatistics.py::marshall or src/m/solve/WriteData.py
#
md.qmu.statistics.nfiles_per_directory = 2
md.qmu.statistics.ndirectories = 5

md.qmu.statistics.method[0]['name'] = 'Histogram'
md.qmu.statistics.method[0]['fields'] = ['Sealevel', 'BslcIce']
md.qmu.statistics.method[0]['steps'] = np.arange(1, 10 + 1).reshape(1, -1)
md.qmu.statistics.method[0]['nbins'] = 20

md.qmu.statistics.addmethod()
md.qmu.statistics.method[1]['name'] = 'MeanVariance'
md.qmu.statistics.method[1]['fields'] = ['Sealevel', 'BslcIce']
md.qmu.statistics.method[1]['steps'] = np.arange(1, 10 + 1).reshape(1, -1)

md.qmu.statistics.addmethod()
md.qmu.statistics.method[2]['name'] = 'SampleSeries'
md.qmu.statistics.method[2]['fields'] = ['Sealevel', 'BslcIce']
md.qmu.statistics.method[2]['steps'] = np.arange(1, 10 + 1).reshape(1, -1)
md.qmu.statistics.method[2]['indices'] = locations.reshape(1, -1)
#}}}

# Run transient Dakota solution
mds = solve(md, 'Transient')

# Run without statistics computations
md.qmu.statistics.method[0]['name'] = 'None'
md = solve(md, 'Transient')

# Compare statistics with our own here
svalues = mds.results.StatisticsSolution[-1].SealevelSamples  # all values at locations

dvalues = np.zeros((md.qmu.method.params.samples, len(locations)))
for i in range(md.qmu.method.params.samples):
    dvalues[i, :] = md.results.dakota.modelresults[i].TransientSolution[-1].Sealevel[locations - 1].flatten()

samplesnorm = np.linalg.norm(dvalues - svalues, 'fro')

# Fields and tolerances to track changes
field_names = ['Samples Norm']
field_tolerances = [1e-13]
field_values = [samplesnorm]
