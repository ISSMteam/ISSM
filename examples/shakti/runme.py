import numpy as np
from model import *
from solve import solve
from hydrologyshakti import hydrologyshakti
from parameterize import parameterize
from export_netCDF import export_netCDF
from transient import transient
from setmask import setmask
from triangle import triangle


steps = [1, 2, 3]

if 1 in steps:
    print('    Step 1: Mesh')
    #Generate unstructured mesh on 1, 000 m square with typical element edge length of 20 m
    md = triangle(model(), './outline.exp', 20)
    export_netCDF(md, 'MoulinMesh.nc')


if 2 in steps:
    print('    Step 2: Parameterization')
    md = loadmodel('MoulinMesh.nc')
    md = setmask(md, '', '')

    # Run parameterization script to set up geometry, velocity, material properties, etc.
    md = parameterize(md, 'moulin.py')

    # HYDROLOGY SPECIFIC PARAMETERIZATION:
    # Change hydrology class to Sommers' SHaKTI model
    md.hydrology = hydrologyshakti()

    # Define initial water head such that water pressure is 50% of ice overburden pressure
    md.hydrology.head = 0.5 * md.materials.rho_ice / md.materials.rho_freshwater * md.geometry.thickness + md.geometry.base

    # Initial subglacial gap height of 0.01m everywhere
    md.hydrology.gap_height = 0.01 * np.ones(md.mesh.numberofelements)

    # Typical bed bump bump spacing (2m)
    md.hydrology.bump_spacing = 2 * np.ones(md.mesh.numberofelements)

    # Typical bed bump height (0.1m)
    md.hydrology.bump_height = 0.1 * np.ones(md.mesh.numberofelements)

    # Define distributed englacial input to the subglacial system (m / yr)
    # Change the value 0.0 to add distributed input
    md.hydrology.englacial_input = 0.0 * np.ones(md.mesh.numberofvertices)

    # Initial Reynolds number (start at Re = 1000 everywhere)
    md.hydrology.reynolds = 1000. * np.ones(md.mesh.numberofelements)

    # Set up atmospheric pressure Type I boundary condition at left edge of
    # domain (outflow, i.e. h = zb at x = xmin)
    md.hydrology.spchead = np.nan * np.ones(md.mesh.numberofvertices)
    pos = np.where(np.logical_and(md.mesh.vertexonboundary, md.mesh.x == np.nanmin(md.mesh.x)))
    md.hydrology.spchead[pos] = md.geometry.base[pos]

    export_netCDF(md, 'MoulinParam.nc')

if 3 in steps:
    print('    Step 3: Solve!')
    md = loadmodel('MoulinParam.nc')

    md.transient = transient.deactivateall(md.transient)
    md.transient.ishydrology = 1

    # Specify that you want to run the model on your current computer
    # Change the number of processors according to your machine (here np = 4)
    md.cluster = generic('np', 2)

    # Define the time stepping scheme: run for 90 days with a time step of 1 hr
    md.timestepping.time_step = 3600. / md.constants.yts  # Time step (in years)
    md.timestepping.final_time = 30. / 365.

    #Add one moulin with steady input at x = 500, y = 500
    DistToCenter = np.sqrt((md.mesh.x - 500)**2 + (md.mesh.y - 500)**2)
    loc = np.where(np.nanmin(DistToCenter) == DistToCenter)
    time = np.arange(0, md.timestepping.final_time + md.timestepping.time_step, md.timestepping.time_step)
    md.hydrology.moulin_input = np.zeros((md.mesh.numberofvertices, len(time)))
    md.hydrology.moulin_input = np.vstack((md.hydrology.moulin_input, time))
    md.hydrology.moulin_input[loc, :] = 4

    # Specify no - flux Type 2 boundary conditions on all edges (except
    # the Type 1 condition set at the outflow above)
    md.hydrology.neumannflux = np.zeros((md.mesh.numberofelements, len(time)))
    md.hydrology.neumannflux = np.vstack((md.hydrology.neumannflux, time))

    md.verbose.solution = 1
    md = solve(md, 'Transient')

    export_netCDF(md, 'MoulinTransient.nc')
