#Test Name: SquareShelfConstrainedSampling
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from solve import *
import numpy as np

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.par')
md = md.sampling.setparameters(md,2e5,1)
md.sampling.seed = 100
md.sampling.phi = np.zeros(md.mesh.numberofvertices, 1)
md.cluster = generic('name', gethostname(), 'np', 1)
md.cluster.np = 1
md = solve(md, 'smp')

#Fields and tolerances to track changes
field_names = ['Sample']
field_tolerances = [1e-13]
field_values = [md.results.SamplingSolution.Sample]
