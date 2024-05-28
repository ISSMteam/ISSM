#Test Name: SquareSheetShelfGroundingLine3dAggressive
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

md = triangle(model(), '../Exp/Square.exp', 350000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.initialization.vx[:] = 0.
md.initialization.vy[:] = 0.
md.geometry.base = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.bed = -700. - abs(md.mesh.y - 500000.) / 1000.
md.geometry.thickness[:] = 1000.
md.geometry.surface = md.geometry.base + md.geometry.thickness
md.smb.mass_balance[:] = 100.
md = setflowequation(md, 'SSA', 'all')
md.transient.isstressbalance = False
md.transient.isgroundingline = True
md.groundingline.migration = 'AggressiveMigration'
md.groundingline.melt_interpolation = 'SubelementMelt1'
md.transient.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation', 'IceMass', 'IceVolumeAboveFloatationScaled', 'GroundedArea', 'FloatingArea', 'GroundedAreaScaled', 'FloatingAreaScaled']
md.mesh.scale_factor = 1.1 * np.ones((md.mesh.numberofvertices))
md.extrude(3, 1.)
md.cluster = generic('name', gethostname(), 'np', 3)

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['Bed1', 'Surface1', 'Thickness1', 'Floatingice1', 'IceVolume1', 'IceVolumeAboveFloatation1', 'IceMass1', 'IceVolumeAboveFloatationScaled1', 'GroundedArea1', 'GroundedAreaScaled1', 'FloatingArea1', 'FloatingAreaScaled1',
               'Bed2', 'Surface2', 'Thickness2', 'Floatingice2', 'IceVolume2', 'IceVolumeAboveFloatation2', 'IceMass2', 'IceVolumeAboveFloatationScaled2', 'GroundedArea2', 'GroundedAreaScaled2', 'FloatingArea2', 'FloatingAreaScaled2',
               'Bed3', 'Surface3', 'Thickness3', 'Floatingice3', 'IceVolume3', 'IceVolumeAboveFloatation3', 'IceMass3', 'IceVolumeAboveFloatationScaled3', 'GroundedArea3', 'GroundedAreaScaled3', 'FloatingArea3', 'FloatingAreaScaled3']

field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-11, 1e-10, 2e-11, 3e-11, 2e-12, 6e-12, 2e-12, 6e-12, 2.5e-12, 2.5e-12, 8e-12, 8e-12,
                    1e-10, 1e-10, 1e-10, 5e-11, 2e-12, 5e-12, 2e-12, 5e-12, 5e-11, 7e-12, 7e-11, 7e-11]

field_values = [md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].MaskOceanLevelset,
                md.results.TransientSolution[0].IceVolume,
                md.results.TransientSolution[0].IceVolumeAboveFloatation,
                md.results.TransientSolution[0].IceMass,
                md.results.TransientSolution[0].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[0].GroundedArea,
                md.results.TransientSolution[0].GroundedAreaScaled,
                md.results.TransientSolution[0].FloatingArea,
                md.results.TransientSolution[0].FloatingAreaScaled,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].MaskOceanLevelset,
                md.results.TransientSolution[1].IceVolume,
                md.results.TransientSolution[1].IceVolumeAboveFloatation,
                md.results.TransientSolution[1].IceMass,
                md.results.TransientSolution[1].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[1].GroundedArea,
                md.results.TransientSolution[1].GroundedAreaScaled,
                md.results.TransientSolution[1].FloatingArea,
                md.results.TransientSolution[1].FloatingAreaScaled,
                md.results.TransientSolution[2].Base,
                md.results.TransientSolution[2].Surface,
                md.results.TransientSolution[2].Thickness,
                md.results.TransientSolution[2].MaskOceanLevelset,
                md.results.TransientSolution[2].IceVolume,
                md.results.TransientSolution[2].IceVolumeAboveFloatation,
                md.results.TransientSolution[2].IceMass,
                md.results.TransientSolution[2].IceVolumeAboveFloatationScaled,
                md.results.TransientSolution[2].GroundedArea,
                md.results.TransientSolution[2].GroundedAreaScaled,
                md.results.TransientSolution[2].FloatingArea,
                md.results.TransientSolution[2].FloatingAreaScaled]
