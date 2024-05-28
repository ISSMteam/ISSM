#Test Name: 79NorthPISMhydro
import numpy as np
from model import *
from setmask import setmask
from triangle import triangle
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from socket import gethostname
from generic import generic


md = triangle(model(), '../Exp/79North.exp', 10000.)
md = setmask(md, '../Exp/79NorthShelf.exp', '')
md = parameterize(md, '../Par/79North.py')
md = setflowequation(md, 'SSA', 'all')

#Hydrology
md.hydrology = hydrologypism()
md.hydrology.drainage_rate = 0.001 * np.ones((md.mesh.numberofvertices))
md.hydrology.watercolumn_max = 20 * np.ones((md.mesh.numberofvertices))
md.initialization.pressure = np.zeros((md.mesh.numberofvertices))
md.initialization.watercolumn = np.zeros((md.mesh.numberofvertices))
md.basalforcings.groundedice_melting_rate = np.arange(0, md.mesh.numberofvertices) + 1
md.transient.ishydrology = 1
md.transient.issmb = 0
md.transient.ismasstransport = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.isgroundingline = 0

md.cluster = generic('name', gethostname(), 'np', 3)
md = solve(md, 'Transient')

#Plot to check result
# plotmodel(md, 'data', md.results.TransientSolution(3).Watercolumn - 3 * (md.materials.rho_freshwater / md.materials.rho_ice * [1:md.mesh.numberofvertices]' - 1))

#Fields and tolerances to track changes
field_names = ['WaterColumn1', 'WaterColumn2', 'WaterColumn3']
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [md.results.TransientSolution[0].Watercolumn,
                md.results.TransientSolution[1].Watercolumn,
                md.results.TransientSolution[2].Watercolumn]
