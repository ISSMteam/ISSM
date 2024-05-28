from model import *
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from generic import generic
from socket import gethostname
from solve import solve

md = model()
md = triangle(md, 'DomainOutline.exp', 100000)
md = setmask(md, 'all', '')
md = parameterize(md, 'Square.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Stressbalance')
