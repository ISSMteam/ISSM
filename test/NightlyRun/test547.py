#Test Name: PigTranFrontalforcingsrignotRunoff
import numpy as np
from calvingvonmises import calvingvonmises
from frontalforcingsrignot import frontalforcingsrignot
from socket import gethostname
from model import *
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from solve import solve
from triangle import triangle

from SMBpddSicopolis import *

md = triangle(model(), '../Exp/Pig.exp', 10000)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md = setflowequation(md, 'SSA', 'all')
md.timestepping.time_step = 1
md.timestepping.final_time = 3

# Separate domain in 2 basins
idbasin = np.zeros((md.mesh.numberofelements,))
iid1 = np.where(md.mesh.x <= -1.6e6)[0]
for ii in range(md.mesh.numberofelements):
    for vertex in range(3):
        if md.mesh.elements[ii][vertex] - 1 in iid1:  # one vertex in basin 1; NOTE: offset because of 1-based vertex indexing
            idbasin[ii] = 1
    if idbasin[ii] == 0:  # no vertex was found in basin 1
        idbasin[ii] = 2

# SMB: PDD scheme providing the runoff that drives the subglacial discharge
# (monthly cycle of -18C to +2C at sea level, 0.4 m/yr precip: ~0.7-1 m/yr of runoff on the low-lying shelves, none at high elevation)
nv = md.mesh.numberofvertices
md.smb = SMBpddSicopolis()
md.smb.initialize(md)
surface = np.maximum(md.geometry.surface, 0)
md.smb.s0p = surface.reshape(-1, 1)
md.smb.s0t = surface.reshape(-1, 1)
temp_ma_present = -18. * np.ones((nv, )) - md.smb.rlaps * surface / 1000.
temp_mj_present = 2. * np.ones((nv, )) - md.smb.rlaps * surface / 1000.
precipitation = 0.4 * np.ones((nv, ))
md.smb.monthlytemperatures = np.zeros((nv, 12))
md.smb.precipitation = np.zeros((nv, 12))
for imonth in range(12):
    md.smb.monthlytemperatures[0:nv, imonth] = md.materials.meltingpoint + temp_ma_present + (temp_mj_present - temp_ma_present) * np.sin((imonth + 1. - 4.) * np.pi / 6.0)
    md.smb.precipitation[0:nv, imonth] = precipitation

# Calving and frontalforcings parameters
md.mask.ice_levelset = 1e4 * (md.mask.ice_levelset + 0.5)
md.calving = calvingvonmises()
md.frontalforcings = frontalforcingsrignot()
md.frontalforcings.num_basins = 2
md.frontalforcings.basin_id = idbasin
md.frontalforcings.isdischargerunoff = 1  # discharge = basin-integrated SMB runoff, subglacial_discharge is not used
md.frontalforcings.thermalforcing = 0.5 * np.ones((nv,))
for elem in range(md.mesh.numberofelements):
    if idbasin[elem] == 2:
        md.frontalforcings.thermalforcing[md.mesh.elements[elem, 0:3] - 1] = 1.5  #NOTE: offset because of 1-based vertex indexing

md.transient.issmb = 1
md.transient.ismovingfront = 1
md.levelset.spclevelset = np.full((nv,), np.nan)
md.levelset.migration_max = 1e10

md.transient.requested_outputs = ['default', 'CalvingMeltingrate', 'SmbRunoff']
md.cluster = generic('name', gethostname(), 'np', 2)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vel1', 'Thickness1', 'MaskIceLevelset1', 'CalvingMeltingrate1', 'SmbRunoff1',
               'Vel2', 'Thickness2', 'MaskIceLevelset2', 'CalvingMeltingrate2', 'SmbRunoff2',
               'Vel3', 'Thickness3', 'MaskIceLevelset3', 'CalvingMeltingrate3', 'SmbRunoff3']
field_tolerances = [
    2e-11, 9e-12, 1e-11, 1e-9, 1e-13,
    2e-11, 2e-11, 1e-11, 1e-9, 1e-13,
    1e-10, 1e-10, 1e-9, 1e-9, 1e-13]
field_values = [
    md.results.TransientSolution[0].Vel,
    md.results.TransientSolution[0].Thickness,
    md.results.TransientSolution[0].MaskIceLevelset,
    md.results.TransientSolution[0].CalvingMeltingrate,
    md.results.TransientSolution[0].SmbRunoff,
    md.results.TransientSolution[1].Vel,
    md.results.TransientSolution[1].Thickness,
    md.results.TransientSolution[1].MaskIceLevelset,
    md.results.TransientSolution[1].CalvingMeltingrate,
    md.results.TransientSolution[1].SmbRunoff,
    md.results.TransientSolution[2].Vel,
    md.results.TransientSolution[2].Thickness,
    md.results.TransientSolution[2].MaskIceLevelset,
    md.results.TransientSolution[2].CalvingMeltingrate,
    md.results.TransientSolution[2].SmbRunoff]
