#Test Name: SquareShelfConstrainedRestartTranSSA2d
from model import *
from socket import gethostname
from triangle import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *
from generic import generic
import copy

md = triangle(model(), '../Exp/Square.exp', 150000.)
md = setmask(md, 'all', '')
md = parameterize(md, '../Par/SquareShelfConstrained.py')
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 1)
md.transient.requested_outputs = ['IceVolume', 'TotalSmb']

md.verbose = verbose('solution', 1)
md.settings.checkpoint_frequency = 4

# time steps and resolution
md.timestepping.final_time = 19
md.settings.output_frequency = 2

md = solve(md, 'Transient')
md2 = copy.deepcopy(md)
md = solve(md, 'Transient', 'restart', 1)

#Fields and tolerances to track changes
field_names = ['Vx1', 'Vy1', 'Vel1', 'TotalSmb1', 'Bed1', 'Surface1', 'Thickness1', 'Volume1', 'Vx2', 'Vy2', 'Vel2', 'TotalSmb2', 'Bed2', 'Surface2', 'Thickness2', 'Volume2', 'Vx3', 'Vy3', 'Vel3', 'TotalSmb3', 'Bed3', 'Surface3', 'Thickness3', 'Volume3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13,
                    1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [md2.results.TransientSolution[6].Vx - md.results.TransientSolution[6].Vx,
                md2.results.TransientSolution[6].Vy - md.results.TransientSolution[6].Vy,
                md2.results.TransientSolution[6].Vel - md.results.TransientSolution[6].Vel,
                md2.results.TransientSolution[6].TotalSmb - md.results.TransientSolution[6].TotalSmb,
                md2.results.TransientSolution[6].Base - md.results.TransientSolution[6].Base,
                md2.results.TransientSolution[6].Surface - md.results.TransientSolution[6].Surface,
                md2.results.TransientSolution[6].Thickness - md.results.TransientSolution[6].Thickness,
                md2.results.TransientSolution[6].IceVolume - md.results.TransientSolution[6].IceVolume,
                md2.results.TransientSolution[7].Vx - md.results.TransientSolution[7].Vx,
                md2.results.TransientSolution[7].Vy - md.results.TransientSolution[7].Vy,
                md2.results.TransientSolution[7].Vel - md.results.TransientSolution[7].Vel,
                md2.results.TransientSolution[7].TotalSmb - md.results.TransientSolution[7].TotalSmb,
                md2.results.TransientSolution[7].Base - md.results.TransientSolution[7].Base,
                md2.results.TransientSolution[7].Surface - md.results.TransientSolution[7].Surface,
                md2.results.TransientSolution[7].Thickness - md.results.TransientSolution[7].Thickness,
                md2.results.TransientSolution[7].IceVolume - md.results.TransientSolution[7].IceVolume,
                md2.results.TransientSolution[8].Vx - md.results.TransientSolution[8].Vx,
                md2.results.TransientSolution[8].Vy - md.results.TransientSolution[8].Vy,
                md2.results.TransientSolution[8].Vel - md.results.TransientSolution[8].Vel,
                md2.results.TransientSolution[8].TotalSmb - md.results.TransientSolution[8].TotalSmb,
                md2.results.TransientSolution[8].Base - md.results.TransientSolution[8].Base,
                md2.results.TransientSolution[8].Surface - md.results.TransientSolution[8].Surface,
                md2.results.TransientSolution[8].Thickness - md.results.TransientSolution[8].Thickness,
                md2.results.TransientSolution[8].IceVolume - md.results.TransientSolution[8].IceVolume]
