#Test Name: SquareShelfConstrainedTranMisfitSurface
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from misfit import *


md = triangle(model(), '../Exp/Square.exp', 180000)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

fake_surface = np.vstack((np.append(np.array(md.geometry.surface) + 100, 1.1),
                          np.append(np.array(md.geometry.surface) + 200, 2.1),
                          np.append(np.array(md.geometry.surface) + 300, 2.5))).T

md.transient.requested_outputs = ['default', 'SurfaceMisfit']
md.outputdefinition.definitions = [misfit(name='SurfaceMisfit',
                                          definitionstring='Outputdefinition1',
                                          model_string='Surface',
                                          observation=fake_surface,
                                          observation_string='SurfaceObservation',
                                          timeinterpolation='nearestneighbor',
                                          weights=np.ones((md.mesh.numberofvertices, 1)),
                                          weights_string='WeightsSurfaceObservation')]

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['SurfaceMisfitFirstStep', 'SurfaceMisfitSecondStep', 'SurfaceMisfitThirdStep']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].SurfaceMisfit,
                md.results.TransientSolution[1].SurfaceMisfit,
                md.results.TransientSolution[2].SurfaceMisfit]
