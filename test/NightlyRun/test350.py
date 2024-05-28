#Test Name: SquareSheetHydrologyShakti
import numpy as np

from frictionshakti import frictionshakti
from socket import gethostname
from model import *
from operator import itemgetter
from parameterize import parameterize
from setflowequation import setflowequation
from setmask import setmask
from solve import solve
from transient import transient
from triangle import triangle

md = triangle(model(), '../Exp/Square.exp', 50000.)
md.mesh.x = md.mesh.x / 1000
md.mesh.y = md.mesh.y / 1000
md = setmask(md, '', '')
md = parameterize(md, '../Par/SquareSheetConstrained.py')
md.transient = transient.deactivateall(md.transient)
md.transient.ishydrology = 1
md = setflowequation(md, 'SSA', 'all')
md.cluster = generic('name', gethostname(), 'np', 2)

#Use hydrology coupled friction law
md.friction = frictionshakti(md.friction)

#Change hydrology class to Shakti' model
md.hydrology = hydrologyshakti()

#Change geometry
md.geometry.base = -.02 * md.mesh.x + 20.
md.geometry.thickness = 300. * np.ones((md.mesh.numberofvertices, ))
md.geometry.bed = md.geometry.base
md.geometry.surface = md.geometry.bed + md.geometry.thickness

#define the initial water head as being such that the water pressure is 50% of the ice overburden pressure
md.hydrology.head = 0.5 * md.materials.rho_ice / md.materials.rho_freshwater * md.geometry.thickness + md.geometry.base
md.hydrology.gap_height = 0.01 * np.ones((md.mesh.numberofelements, ))
md.hydrology.bump_spacing = 2 * np.ones((md.mesh.numberofelements, ))
md.hydrology.bump_height = 0.05 * np.ones((md.mesh.numberofelements, ))
md.hydrology.englacial_input = 0.5 * np.ones((md.mesh.numberofvertices, ))
md.hydrology.reynolds = 1000. * np.ones((md.mesh.numberofelements, ))
md.hydrology.spchead = float('NaN') * np.ones((md.mesh.numberofvertices, ))
pos = np.intersect1d(np.array(np.where(md.mesh.vertexonboundary)), np.array(np.where(md.mesh.x == 1000)))
md.hydrology.spchead[pos] = md.geometry.base[pos]

#Define velocity
md.initialization.vx = 1e-6 * md.constants.yts * np.ones((md.mesh.numberofvertices, ))
md.initialization.vy = np.zeros((md.mesh.numberofvertices, ))

md.timestepping.time_step = 3. * 3600. / md.constants.yts
md.timestepping.final_time = .5 / 365.
md.materials.rheology_B = (5e-25)**(-1. / 3.) * np.ones((md.mesh.numberofvertices, ))

#Add one moulin and Neumann BC, varying in time
a = np.sqrt((md.mesh.x - 500.)**2 + (md.mesh.y - 500.)**2)
pos = min(enumerate(a), key=itemgetter(1))[0]
time = np.arange(0, md.timestepping.final_time + 1, md.timestepping.time_step)
md.hydrology.moulin_input = np.zeros((md.mesh.numberofvertices + 1, np.size(time)))
md.hydrology.moulin_input[-1, :] = time
md.hydrology.moulin_input[pos, :] = 5. * (1. - np.sin(2. * np.pi / (1. / 365.) * time))
md.hydrology.neumannflux = np.zeros((md.mesh.numberofelements + 1, np.size(time)))
md.hydrology.neumannflux[-1, :] = time
segx = md.mesh.x[md.mesh.segments[:, 0] - 1]
segy = md.mesh.y[md.mesh.segments[:, 0] - 1]
posA = np.intersect1d(np.intersect1d(np.array(np.where(segx < 1.)), np.array(np.where(segy > 400.))), np.array(np.where(segy < 600.)))
pos = (md.mesh.segments[posA] - 1)[:, 2]
md.hydrology.neumannflux[pos, :] = np.tile(0.05 * (1. - np.sin(2. * np.pi / (1. / 365.) * time)), (len(pos), 1))

md = solve(md, 'Transient')

#Fields and tolerances to track changes
field_names = ['HydrologyHead1', 'HydrologyGapHeight1',
               'HydrologyHead2', 'HydrologyGapHeight2',
               'HydrologyHead3', 'HydrologyGapHeight3',
               'HydrologyHead4', 'HydrologyGapHeight4']
field_tolerances = [1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-13,
                    1e-13, 1e-12]
field_values = [md.results.TransientSolution[0].HydrologyHead,
                md.results.TransientSolution[0].HydrologyGapHeight,
                md.results.TransientSolution[1].HydrologyHead,
                md.results.TransientSolution[1].HydrologyGapHeight,
                md.results.TransientSolution[2].HydrologyHead,
                md.results.TransientSolution[2].HydrologyGapHeight,
                md.results.TransientSolution[3].HydrologyHead,
                md.results.TransientSolution[3].HydrologyGapHeight]
