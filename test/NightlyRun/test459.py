#Test Name: SquareSheetShelfHydrologyGlaDSTransitionModel
import numpy as np
from model import *
from triangle import triangle
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from SetIceSheetBC import SetIceSheetBC
from generic import generic
from cuffey import cuffey

# Create model
md = triangle(model(), '../Exp/Square.exp', 50000.)
md = setmask(md, '../Exp/SquareShelf.exp', '')
md = parameterize(md, '../Par/SquareSheetShelf.py')
md.mesh.x = md.mesh.x / 100
md.mesh.y = md.mesh.y / 100
md.miscellaneous.name = 'testChannels'

# Miscellaneous
md = setflowequation(md, 'SSA', 'all')

# Some constants
md.constants.g = 9.8
md.materials.rho_ice = 910

# Define initial conditions
md.initialization.vx = 1.0e-6 * md.constants.yts * np.ones((md.mesh.numberofvertices))
md.initialization.vy = np.zeros((md.mesh.numberofvertices))
md.initialization.temperature = (273. - 20.) * np.ones((md.mesh.numberofvertices))
md.initialization.watercolumn = 0.03 * np.ones((md.mesh.numberofvertices))
md.initialization.hydraulic_potential = md.materials.rho_ice * md.constants.g * md.geometry.thickness

# Materials
md.materials.rheology_B = (5e-25)**(-1. / 3.) * np.ones((md.mesh.numberofvertices))
md.materials.rheology_n = 3. * np.ones((md.mesh.numberofelements))

# Friction
md.friction.coefficient = np.zeros((md.mesh.numberofvertices))
md.friction.p = np.ones((md.mesh.numberofelements))
md.friction.q = np.ones((md.mesh.numberofelements))

# Boundary conditions:
md = SetIceSheetBC(md)

md.transient = transient.deactivateall(md.transient)
md.transient.ishydrology = 1
md.transient.isgroundingline = 1

# Set numerical conditions
md.timestepping.time_step = 0.1 / 365
md.timestepping.final_time = 0.4 / 365

# Change hydrology class to Glads model
md.hydrology = hydrologyglads()
md.hydrology.maxiter = 10 # Make sure it runs quickly...
md.hydrology.ischannels = 1
md.hydrology.istransition = 1
md.hydrology.omega = 1 / 2000.
md.hydrology.sheet_alpha = 3. / 2.
md.hydrology.sheet_beta = 3. / 2.
md.hydrology.englacial_void_ratio = 1.e-5
md.hydrology.moulin_input = np.zeros((md.mesh.numberofvertices))
md.hydrology.neumannflux = np.zeros((md.mesh.numberofelements))
md.hydrology.bump_height = 1.e-1 * np.ones((md.mesh.numberofvertices))
md.hydrology.sheet_conductivity = 1.e-2 * np.ones((md.mesh.numberofvertices))
md.hydrology.channel_conductivity = 5.e-2 * np.ones((md.mesh.numberofvertices))
md.hydrology.rheology_B_base = cuffey(273.15 - 2) * np.ones((md.mesh.numberofvertices))

# BCs for hydrology
pos = np.nonzero(np.logical_and(md.mesh.y <= 100, md.mesh.vertexonboundary))
md.hydrology.spcphi = np.nan * np.ones((md.mesh.numberofvertices))
md.hydrology.spcphi[pos] = md.materials.rho_ice * md.constants.g * md.geometry.thickness[pos]

md.cluster = generic('np', 2)
md.hydrology.requested_outputs = ['default', 'TotalHydrologyBasalFlux']
md = solve(md, 'Transient') # Or 'tr'

# Fields and tolerances to track changes
field_names = ['HydrologySheetThickness1', 'HydraulicPotential1', 'ChannelArea1', 'TotalHydrologyBasalFlux',
               'HydrologySheetThickness2', 'HydraulicPotential2', 'ChannelArea2', 'TotalHydrologyBasalFlux',
               'HydrologySheetThickness3', 'HydraulicPotential3', 'ChannelArea3', 'TotalHydrologyBasalFlux',
               'HydrologySheetThickness4', 'HydraulicPotential4', 'ChannelArea4', 'TotalHydrologyBasalFlux']
field_tolerances = [1e-14, 8e-14, 3e-12, 1e-13,
                    1e-14, 8e-14, 3e-12, 1e-13,
                    1e-14, 8e-14, 3e-12, 1e-13,
                    1e-14, 9e-14, 3e-12, 1e-13]
field_values = [md.results.TransientSolution[0].HydrologySheetThickness,
                md.results.TransientSolution[0].HydraulicPotential,
                md.results.TransientSolution[0].ChannelArea,
                md.results.TransientSolution[0].TotalHydrologyBasalFlux,
                md.results.TransientSolution[1].HydrologySheetThickness,
                md.results.TransientSolution[1].HydraulicPotential,
                md.results.TransientSolution[1].ChannelArea,
                md.results.TransientSolution[1].TotalHydrologyBasalFlux,
                md.results.TransientSolution[2].HydrologySheetThickness,
                md.results.TransientSolution[2].HydraulicPotential,
                md.results.TransientSolution[2].ChannelArea,
                md.results.TransientSolution[2].TotalHydrologyBasalFlux,
                md.results.TransientSolution[3].HydrologySheetThickness,
                md.results.TransientSolution[3].HydraulicPotential,
                md.results.TransientSolution[3].ChannelArea,
                md.results.TransientSolution[3].TotalHydrologyBasalFlux]
