# Mismip3D experiment with AMR using BAMG
from bamg import bamg
from model import *
from export_netCDF import export_netCDF
from setmask import setmask
from loadmodel import loadmodel
from parameterize import parameterize
from generic import generic
from exportVTK import exportVTK
from setflowequation import setflowequation
from solve import solve
from plotmodel import plotmodel


steps = [4]

if 1 in steps:
    print('   Step 1: Coarse mesh')

    #Generate an unstructured coarse mesh on the MISMIP domain with typical element edge length equal to 10, 000 m
    md = bamg(model(), 'domain', './domain.exp', 'hmax', 10000, 'splitcorners', 1)

    export_netCDF(md, 'AMRCoarseMesh.nc')

if 2 in steps:
    print('   Step 2: Parameterization')

    md = loadmodel('AMRCoarseMesh.nc')

    md = setmask(md, '', '')

    # Run parameterization script to set up geometry, inital velocity, material properties, etc.
    md = parameterize(md, './mismip.py')

    # Set the AMR properties and the refinement criteria
    # Here, we are refining around the grounding line
    # We impose the element resolution at grounding equal to 1000 m (1 km)
    # The criterion used is the element distance to the grounding line
    # The distance used here is 10000 m (10 km), used in both side around the grouding line (upstream and downstream)
    md.amr.groundingline_resolution = 1000
    md.amr.groundingline_distance = 10000
    md.amr.hmin = 1000  # the same resolution used around the grounding line
    md.amr.hmax = 10000  # the same coase resolution used to generate the coarse mesh
    md.amr.gradation = 1.7  # this controls the ratio between two consecutive edges
    md.amr.fieldname = 'None'  # no field used here
    md.amr.keepmetric = 0  # no field, no metric

    export_netCDF(md, 'AMRParam.nc')

if 3 in steps:
    print('   Step 3: Solve!')

    md = loadmodel('AMRParam.nc')

    # Run transient with adaptive mesh refinement
    md.timestepping.time_step = 1
    md.timestepping.final_time = 500   # here, as example, only 500 yr.
    md.settings.output_frequency = 10  # here, save results every 10 yr
    md.stressbalance.maxiter = 30
    md.stressbalance.abstol = np.nan
    md.stressbalance.restol = 1
    md.settings.solver_residue_threshold = 1e-2  # relaxing (the first stress balance solver iteration presents values higher than the original threshold. This probably happens because the initial velocity is set to one).
    md.verbose = verbose('convergence', False, 'solution', True)

    # Specify that you want to run the model on your current (local host) computer
    # Change the number of processors according to your machine (here np = 2)
    md.cluster = generic('np', 2)

    # Set the AMR frequency, i.e., can be 1 or larger depending on how often the mesh needs to be updated
    md.transient.amr_frequency = 1  # here, we are refining the mesh in every time step

    # Set the flow equation (SSA) and run
    md = setflowequation(md, 'SSA', 'all')
    md = solve(md, 'Transient')

    # Print the solutions and the mesh in VTK format (needs ParaView:    https: / / www.paraview.org)
    exportVTK('./VTKpy', md)
    export_netCDF(md, 'AMRTransient.nc')
