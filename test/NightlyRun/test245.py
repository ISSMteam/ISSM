#Test Name: SquareShelfTranIspddSicopolisSSA2d
import numpy as np
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

from SMBpddSicopolis import *

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelf.py')

# Use of SMBpddSicopolis
md.smb = SMBpddSicopolis()
# initalize pdd fields
md.smb.initialize(md)
md.smb.s0p = md.geometry.surface.reshape(-1, 1)
md.smb.s0t = md.geometry.surface.reshape(-1, 1)


md.smb.monthlytemperatures = np.empty((md.mesh.numberofvertices, 12))
md.smb.precipitation = np.empty((md.mesh.numberofvertices, 12))
temp_ma_present = -10. * np.ones((md.mesh.numberofvertices, )) - md.smb.rlaps * md.geometry.surface / 1000.
temp_mj_present = 10. * np.ones((md.mesh.numberofvertices, )) - md.smb.rlaps * md.geometry.surface / 1000.
precipitation = 5. * np.ones((md.mesh.numberofvertices, ))
for imonth in range(12):
    md.smb.monthlytemperatures[0:md.mesh.numberofvertices, imonth] = md.materials.meltingpoint + temp_ma_present + (temp_mj_present - temp_ma_present) * np.sin((imonth + 1. - 4.) * np.pi / 6.0)
    md.smb.precipitation[0:md.mesh.numberofvertices, imonth] = precipitation

# time steps and resolution
md.timestepping.time_step = 1
md.settings.output_frequency = 1
md.timestepping.final_time = 2

md.transient.issmb = 1
md.transient.ismasstransport = 1
md.transient.isstressbalance = 0
md.transient.isthermal = 0

md.transient.requested_outputs = ['default', 'TemperaturePDD']
md.cluster = generic('name', gethostname(), 'np', 1)  # 3 for the cluster
md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['TemperaturePDD1', 'SmbMassBalance1', 'TemperaturePDD2', 'SmbMassBalance2']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].TemperaturePDD,
                md.results.TransientSolution[0].SmbMassBalance,
                md.results.TransientSolution[1].TemperaturePDD,
                md.results.TransientSolution[1].SmbMassBalance]
