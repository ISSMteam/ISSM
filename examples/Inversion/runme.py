import numpy as np
from model import *
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from generic import generic
from solve import solve
from plotmodel import plotmodel
from export_netCDF import export_netCDF
from m1qn3inversion import m1qn3inversion
from verbose import verbose
from loadmodel import loadmodel
from cuffey import cuffey
steps = [4]
Clims = [1.3 * 1e8, 1.9 * 1e8]

if 1 in steps:
    #Generate observations
    md = model()
    md = triangle(md, 'DomainOutline.exp', 100000)
    md = setmask(md, 'all', '')
    md = parameterize(md, 'Square.py')
    md = setflowequation(md, 'SSA', 'all')
    md.cluster = generic('np', 2)
    md = solve(md, 'Stressbalance')
    plotmodel(md, 'axis#all', 'tight', 'data', md.materials.rheology_B, 'caxis', Clims, 'title', '"True" B',
              'data', md.results.StressbalanceSolution.Vel, 'title', '"observed velocities"')
    export_netCDF(md, 'model1.nc')

if 2 in steps:
    #Modify rheology, now constant
    md = loadmodel('model1.nc')
    md.materials.rheology_B[:] = 1.8 * 1e8

    #results of previous run are taken as observations
    md.inversion = m1qn3inversion()
    md.inversion.vx_obs = md.results.StressbalanceSolution.Vx
    md.inversion.vy_obs = md.results.StressbalanceSolution.Vy
    md.inversion.vel_obs = md.results.StressbalanceSolution.Vel

    md = solve(md, 'Stressbalance')
    plotmodel(md, 'axis#all', 'tight', 'data', md.materials.rheology_B, 'caxis', Clims, 'title', 'B first guess',
              'data', md.results.StressbalanceSolution.Vel, 'title', 'modeled velocities')
    export_netCDF(md, 'model2.nc')


if 3 in steps:
    #invert for ice rigidity
    md = loadmodel('model2.nc')

    #Set up inversion parameters
    maxsteps = 20
    md.inversion.iscontrol = 1
    md.inversion.control_parameters = ['MaterialsRheologyBbar']
    md.inversion.maxsteps = maxsteps
    md.inversion.cost_functions = [101]
    md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 1))
    md.inversion.min_parameters = cuffey(273) * np.ones((md.mesh.numberofvertices, 1))
    md.inversion.max_parameters = cuffey(200) * np.ones((md.mesh.numberofvertices, 1))

    #Go solve!
    md.verbose = verbose(0)
    md = solve(md, 'Stressbalance')
    plotmodel(md, 'axis#all', 'tight', 'data', md.results.StressbalanceSolution.MaterialsRheologyBbar, 'caxis', Clims, 'title', 'inferred B',
              'data', md.results.StressbalanceSolution.Vel, 'title', 'modeled velocities')


if 4 in steps:
    #invert for ice rigidity
    md = loadmodel('model2.nc')

    #Set up inversion parameters
    maxsteps = 20
    md.inversion.iscontrol = 1
    md.inversion.control_parameters = ['MaterialsRheologyBbar']
    md.inversion.maxsteps = maxsteps
    md.inversion.cost_functions = [101, 502]
    md.inversion.cost_functions_coefficients = np.ones((md.mesh.numberofvertices, 2))
    md.inversion.cost_functions_coefficients[:, 1] = 1e-16
    md.inversion.min_parameters = cuffey(273) * np.ones((md.mesh.numberofvertices, 1))
    md.inversion.max_parameters = cuffey(200) * np.ones((md.mesh.numberofvertices, 1))

    #Go solve!
    md.verbose = verbose(0)
    md = solve(md, 'Stressbalance')
    plotmodel(md, 'axis#all', 'tight', 'data', md.results.StressbalanceSolution.MaterialsRheologyBbar, 'caxis', Clims, 'title', 'inferred B',
              'data', md.results.StressbalanceSolution.Vel, 'title', 'modeled velocities')
