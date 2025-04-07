#Test Name: SquareSheetHydrologyGlaDSSheetInPhi
import numpy as np
from model import *
from triangle import triangle
from setmask import setmask
from setflowequation import setflowequation
from parameterize import parameterize
from solve import solve
from SetIceSheetBC import SetIceSheetBC
from generic import generic

# Create model
md = triangle(model(), '../Exp/Square.exp', 50000.)
md.mesh.x = md.mesh.x / 100
md.mesh.y = md.mesh.y / 100
md.miscellaneous.name = 'testChannels'

# Miscellaneous
md = setmask(md, '', '') # Everywhere grounded
md = setflowequation(md, 'SSA', 'all')
md.stressbalance.maxiter = 10 # Make sure it runs quickly...

# Some constants
md.constants.g = 9.8
md.materials.rho_ice = 910

# Geometry
md.geometry.surface = -0.02 * md.mesh.x + 320
md.geometry.bed = np.zeros((md.mesh.numberofvertices))
md.geometry.base = md.geometry.bed
md.geometry.thickness = md.geometry.surface - md.geometry.bed

# Define initial conditions
md.initialization.vx = 1.0e-6 * md.constants.yts * np.ones((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.initialization.watercolumn = 0.03 * np.ones((md.mesh.numberofvertices))
md.initialization.hydraulic_potential = md.materials.rho_ice * md.constants.g * md.geometry.thickness

#cMaterials
md.materials.rheology_B = (5e-25)**(-1./3.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

#cFriction
md.friction.coefficient = np.zeros((md.mesh.numberofvertices))
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))
#md.friction.coupling = 0

#Bcoundary conditions:
md = SetIceSheetBC(md)

md.inversion.iscontrol = 0
md.transient = transient.deactivateall(md.transient)
md.transient.ishydrology = 1

# Set numerical conditions
md.timestepping.time_step = 0.1 / 365
md.timestepping.final_time = 0.4 / 365

#Change hydrology class to Glads model
md.hydrology = hydrologyglads()
md.hydrology.ischannels = 1
md.hydrology.isincludesheetthickness = 1
md.hydrology.englacial_void_ratio = 1.e-5
md.hydrology.moulin_input = np.zeros((md.mesh.numberofvertices))
md.hydrology.neumannflux = np.zeros((md.mesh.numberofelements))
md.hydrology.bump_height = 1.e-1 * np.ones((md.mesh.numberofvertices))
md.hydrology.sheet_conductivity = 1.e-3 * np.ones((md.mesh.numberofvertices))
md.hydrology.channel_conductivity = 5.e-2 * np.ones((md.mesh.numberofvertices))
md.hydrology.rheology_B_base = 8.378836055370960e+07 * np.ones((md.mesh.numberofvertices)) 

# BCs for hydrology
pos = np.where(np.logical_and(md.mesh.x == 100, md.mesh.vertexonboundary))
md.hydrology.spcphi = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcphi[pos] = md.materials.rho_ice * md.constants.g * md.geometry.thickness[pos]

md.cluster = generic('np', 2)
md = solve(md, 'Transient') # Or 'tr'

# Fields and tolerances to track changes
field_names = ['HydrologySheetThickness1', 'HydraulicPotential1', 'ChannelArea1',
               'HydrologySheetThickness2', 'HydraulicPotential2', 'ChannelArea2',
               'HydrologySheetThickness3', 'HydraulicPotential3', 'ChannelArea3',
               'HydrologySheetThickness4', 'HydraulicPotential4', 'ChannelArea4']
field_tolerances = [1e-14, 7e-14, 3e-12,
                    1e-14, 7e-14, 3e-12,
                    1e-14, 7e-14, 3e-12,
                    1e-14, 8e-14, 3e-12]
field_values = [md.results.TransientSolution[0].HydrologySheetThickness,
                md.results.TransientSolution[0].HydraulicPotential,
                md.results.TransientSolution[0].ChannelArea,
                md.results.TransientSolution[1].HydrologySheetThickness,
                md.results.TransientSolution[1].HydraulicPotential,
                md.results.TransientSolution[1].ChannelArea,
                md.results.TransientSolution[2].HydrologySheetThickness,
                md.results.TransientSolution[2].HydraulicPotential,
                md.results.TransientSolution[2].ChannelArea,
                md.results.TransientSolution[3].HydrologySheetThickness,
                md.results.TransientSolution[3].HydraulicPotential,
                md.results.TransientSolution[3].ChannelArea]
