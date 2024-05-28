#Test Name: SquareShelfConstrainedTherTranAdolcReverseVsForward
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from independent import *
from dependent import *
from solve import *
from issmgslsolver import issmgslsolver


#test reverse scalar vs forward vectorial drivers in ADOLC, using the test3009 setup, equivalent to test109 setup.
md = triangle(model(), '../Exp/Square.exp', 100000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md.toolkits.DefaultAnalysis = issmgslsolver()

md.autodiff.isautodiff = True
md.verbose.autodiff = True

#first run scalar reverse mode:
md.autodiff.independents = [independent('name', 'md.geometry.thickness', 'type', 'vertex', 'nods', md.mesh.numberofvertices)]
md.autodiff.dependents = [dependent('name', 'MaxVel', 'type', 'scalar', 'fos_reverse_index', 1)]
md.autodiff.driver = 'fos_reverse'

md = solve(md, 'Transient')

#recover jacobian:
jac_reverse = md.results.TransientSolution[0].AutodiffJacobian

#now run vectorial forward mode
md.autodiff.independents = [independent('name', 'md.geometry.thickness', 'type', 'vertex', 'nods', md.mesh.numberofvertices, 'fov_forward_indices', np.arange(1, md.mesh.numberofvertices + 1))]
md.autodiff.dependents = [dependent('name', 'MaxVel', 'type', 'scalar')]
md.autodiff.driver = 'fov_forward'

md = solve(md, 'Transient')

#recover jacobian:
jac_forward = md.results.TransientSolution[0].AutodiffJacobian

#Fields and tolerances to track changes
field_names = ['Jac Forward', 'Jac Reverse', 'Jac Forward - Reverse']
field_tolerances = [1e-8, 1e-8, 5e-6]
field_values = [jac_forward, jac_reverse, jac_forward - jac_reverse]
