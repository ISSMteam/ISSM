#Test Name: SquareSheetShelfDiadSSA3dDakotaAreaAverage
# this test may crash
import numpy as np

from AreaAverageOntoPartition import *
from model import *
from parameterize import *
from partitioner import *
from setflowequation import *
from setmask import *
from socket import gethostname
from solve import *
from triangle import *


#test partitioning, and partition averaging
# python cannot handle resolutions greater than 30010
md = triangle(model(), '../Exp/Square.exp', 30000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 3)

#partitioning
npart = 100

# Partitioner seemed to generate the following message,
#
#	corrupted size vs. prev_size
#	Aborted (core dumped)
#
# In testing binaries,
#
#	malloc(): invalid next size (unsorted)
#
# TODO:
# - Run valgrind and fix the above
#
partition = partitioner(md, 'package', 'chaco', 'npart', npart)[0] - 1

vector = np.arange(1, 1 + md.mesh.numberofvertices, 1).reshape(-1, 1)
vector_on_partition = AreaAverageOntoPartition(md, vector, partition)
vector_on_nodes = vector_on_partition[partition]

field_names = ['vector_on_nodes']
field_tolerances = [1e-11]
field_values = [vector_on_nodes]
