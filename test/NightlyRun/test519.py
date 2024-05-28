#Test Name: PigTranMOLHO2d
import numpy as np
from model import *
from socket import gethostname
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from SetMOLHOBC import SetMOLHOBC
from generic import generic

md = triangle(model(), '../Exp/Pig.exp', 20000.)
md = setmask(md, '../Exp/PigShelves.exp', '../Exp/PigIslands.exp')
md = parameterize(md, '../Par/Pig.py')
md = setflowequation(md, 'MOLHO', 'all')
md.mesh.scale_factor = 0.9 * np.ones((md.mesh.numberofvertices))
md.transient.requested_outputs = ['default', 'IceVolume', 'IceVolumeScaled', 'GroundedArea', 'GroundedAreaScaled', 'FloatingArea', 'FloatingAreaScaled', 'TotalSmb', 'TotalSmbScaled', 'TotalFloatingBmb', 'TotalFloatingBmbScaled']
md.cluster = generic('name', gethostname(), 'np', 3)
md = SetMOLHOBC(md)
md = solve(md, 'Transient')

# Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'Pressure1',
               'Bed1', 'Surface1', 'Thickness1',
               'IceVolume1', 'IceVolumeScaled1',
               'GroundedArea1', 'GroundedAreaScaled1',
               'FloatingArea1', 'FloatingAreaScaled1',
               'TotalSmb1', 'TotalSmbScaled1',
               'TotalFloatingBmb1', 'TotalFloatingBmbScaled1',
               'Vx2', 'Vy2', 'Vel2', 'Pressure2',
               'Bed2', 'Surface2', 'Thickness2',
               'IceVolume2', 'IceVolumeScaled2',
               'GroundedArea2', 'GroundedAreaScaled2',
               'FloatingArea2', 'FloatingAreaScaled2',
               'TotalSmb2', 'TotalSmbScaled2',
               'TotalFloatingBmb2', 'TotalFloatingBmbScaled2']
field_tolerances = [1e-12, 2e-12, 2e-12, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-12, 1e-12,
                    1e-12, 1e-13,
                    1e-13, 1e-13,
                    4e-13, 4e-13, 4e-13, 1e-13,
                    1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13]
field_values = [md.results.TransientSolution[0].Vx,
                md.results.TransientSolution[0].Vy,
                md.results.TransientSolution[0].Vel,
                md.results.TransientSolution[0].Pressure,
                md.results.TransientSolution[0].Base,
                md.results.TransientSolution[0].Surface,
                md.results.TransientSolution[0].Thickness,
                md.results.TransientSolution[0].IceVolume,
                md.results.TransientSolution[0].IceVolumeScaled,
                md.results.TransientSolution[0].GroundedArea,
                md.results.TransientSolution[0].GroundedAreaScaled,
                md.results.TransientSolution[0].FloatingArea,
                md.results.TransientSolution[0].FloatingAreaScaled,
                md.results.TransientSolution[0].TotalSmb,
                md.results.TransientSolution[0].TotalSmbScaled,
                md.results.TransientSolution[0].TotalFloatingBmb,
                md.results.TransientSolution[0].TotalFloatingBmbScaled,
                md.results.TransientSolution[1].Vx,
                md.results.TransientSolution[1].Vy,
                md.results.TransientSolution[1].Vel,
                md.results.TransientSolution[1].Pressure,
                md.results.TransientSolution[1].Base,
                md.results.TransientSolution[1].Surface,
                md.results.TransientSolution[1].Thickness,
                md.results.TransientSolution[1].IceVolume,
                md.results.TransientSolution[1].IceVolumeScaled,
                md.results.TransientSolution[1].GroundedArea,
                md.results.TransientSolution[1].GroundedAreaScaled,
                md.results.TransientSolution[1].FloatingArea,
                md.results.TransientSolution[1].FloatingAreaScaled,
                md.results.TransientSolution[1].TotalSmb,
                md.results.TransientSolution[1].TotalSmbScaled,
                md.results.TransientSolution[1].TotalFloatingBmb,
                md.results.TransientSolution[1].TotalFloatingBmbScaled]
