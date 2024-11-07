# Test Name: FourierLoveKernels
# Homogenous Earth, for which analytic solutions exist.
# Love kernels for degree 2, 20, 200 (tested against analytic solns).
# Skip benchmarking for the inner-most interface.


from socket import gethostname

import numpy as np

from generic import generic
from materials import *
from model import *
from solve import *
from yi_analytic_homogenous import love_analytic as yi_analytic_homogenous

# Set validation = 1 for comparing against the analytic solutions
validation = 0

# For volumetric potential
md = model()

md.materials = materials('litho')
cst = 365.25 * 24 * 3600 * 1000

md.materials.numlayers = 40
md.love.forcing_type = 9

md.materials.density = np.zeros((md.materials.numlayers, 1)) + 5511
md.materials.lame_mu = np.zeros((md.materials.numlayers, 1)) + 0.75e11
md.materials.lame_lambda = np.zeros((md.materials.numlayers, 1)) + 5e17
md.materials.issolid = np.ones((md.materials.numlayers, 1))
md.materials.rheologymodel = np.zeros((md.materials.numlayers, 1))

# The following isn't used here but needs to have arrays of consistent size with the rest of the materials
md.materials.viscosity = np.zeros((md.materials.numlayers, 1)) + 1e21

md.materials.radius = np.linspace(10e3, 6371e3, md.materials.numlayers + 1).reshape(-1, 1)
md.love.g0 = 9.8134357285509388 # directly grabbed from fourierlovesolver for this particular case

md.love.allow_layer_deletion = 0
md.love.frequencies = 0
md.love.nfreq = 1
md.love.istemporal = 0

md.love.sh_nmin = 2
md.love.sh_nmax = 200
md.love.love_kernels = 1

md.miscellaneous.name = 'kernels'
md.cluster = generic('name', gethostname(), 'np', 3)
md.verbose = verbose('111111101')

md = solve(md,'lv')

# Save yi's for all layers except for the inner-most one, at select degrees.
degrees = [1, 19, 199] # we archive solutions for degrees 2, 20, 200
kernels = np.reshape(md.results.LoveSolution.LoveKernels, (md.love.sh_nmax + 1, md.love.nfreq, md.materials.numlayers + 1, 6), order='F')

# Extract love kernels #{{{
# degree 2
y1_tidal_degree002 = kernels[degrees[0]+1,0,1:,0].squeeze()
y2_tidal_degree002 = kernels[degrees[0]+1,0,1:,1].squeeze()
y3_tidal_degree002 = kernels[degrees[0]+1,0,1:,2].squeeze()
y4_tidal_degree002 = kernels[degrees[0]+1,0,1:,3].squeeze()
y5_tidal_degree002 = kernels[degrees[0]+1,0,1:,4].squeeze()
y6_tidal_degree002 = kernels[degrees[0]+1,0,1:,5].squeeze()

# degree 20
y1_tidal_degree020 = kernels[degrees[1]+1,0,1:,0].squeeze()
y2_tidal_degree020 = kernels[degrees[1]+1,0,1:,1].squeeze()
y3_tidal_degree020 = kernels[degrees[1]+1,0,1:,2].squeeze()
y4_tidal_degree020 = kernels[degrees[1]+1,0,1:,3].squeeze()
y5_tidal_degree020 = kernels[degrees[1]+1,0,1:,4].squeeze()
y6_tidal_degree020 = kernels[degrees[1]+1,0,1:,5].squeeze()

# degree 200
y1_tidal_degree200 = kernels[degrees[2]+1,0,1:,0].squeeze()
y2_tidal_degree200 = kernels[degrees[2]+1,0,1:,1].squeeze()
y3_tidal_degree200 = kernels[degrees[2]+1,0,1:,2].squeeze()
y4_tidal_degree200 = kernels[degrees[2]+1,0,1:,3].squeeze()
y5_tidal_degree200 = kernels[degrees[2]+1,0,1:,4].squeeze()
y6_tidal_degree200 = kernels[degrees[2]+1,0,1:,5].squeeze()
#}}}

# Validate tidal potential solutions against the analytic solutions
if validation:
    param = {}
    param['rho'] = md.materials.density[0]
    param['mu'] = md.materials.lame_mu[1]
    param['G'] = 6.67e-11
    param['radius'] = md.materials.radius
    param['g0'] = md.love.g0
    param['source'] = md.love.forcing_type

    # Check against analytic solutions at the following degrees
    degrees_analytic = [2, 4, 6, 16, 32]
    for jj in range(0, len(degrees_analytic)):
        param['degree'] = degrees_analytic[jj]
        y_temp = yi_analytic_homogenous(param)
        if 'y_ana' not in locals(): # NOTE: Probably not a very Pythonic thing to do, but translated from MATLAB
            y_ana = np.zeros((len(degrees_analytic), np.shape(y_temp)[0], np.shape(y_temp)[1]))
        y_ana[jj,:,:] = y_temp

    y1_ana = np.squeeze(y_ana[:,:,0])
    y2_ana = np.squeeze(y_ana[:,:,1])
    y3_ana = np.squeeze(y_ana[:,:,2])
    y4_ana = np.squeeze(y_ana[:,:,3])
    y5_ana = np.squeeze(y_ana[:,:,4])
    y6_ana = np.squeeze(y_ana[:,:,5])

    depth = (np.max(param['radius']) - param['radius']) / 1000 # km

    kernels = np.reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax + 1, md.love.nfreq, md.materials.numlayers + 1, 6], order='F')

    y1 = np.squeeze(kernels[:,0,:,0])
    y2 = np.squeeze(kernels[:,0,:,1])
    y3 = np.squeeze(kernels[:,0,:,2])
    y4 = np.squeeze(kernels[:,0,:,3])
    y5 = np.squeeze(kernels[:,0,:,4])
    y6 = np.squeeze(kernels[:,0,:,5])

    # TODO: Generate plots similar to those created in MATLAB test

# For surface load
md.love.forcing_type = 11
md = solve(md,'lv')
kernels = np.reshape(md.results.LoveSolution.LoveKernels, [md.love.sh_nmax + 1, md.love.nfreq, md.materials.numlayers + 1, 6], order='F')

# Extract love kernels #{{{
# degree 2
y1_loading_degree002 = kernels[degrees[0]+1,0,1:,0].squeeze()
y2_loading_degree002 = kernels[degrees[0]+1,0,1:,1].squeeze()
y3_loading_degree002 = kernels[degrees[0]+1,0,1:,2].squeeze()
y4_loading_degree002 = kernels[degrees[0]+1,0,1:,3].squeeze()
y5_loading_degree002 = kernels[degrees[0]+1,0,1:,4].squeeze()
y6_loading_degree002 = kernels[degrees[0]+1,0,1:,5].squeeze()

# degree 20
y1_loading_degree020 = kernels[degrees[1]+1,0,1:,0].squeeze()
y2_loading_degree020 = kernels[degrees[1]+1,0,1:,1].squeeze()
y3_loading_degree020 = kernels[degrees[1]+1,0,1:,2].squeeze()
y4_loading_degree020 = kernels[degrees[1]+1,0,1:,3].squeeze()
y5_loading_degree020 = kernels[degrees[1]+1,0,1:,4].squeeze()
y6_loading_degree020 = kernels[degrees[1]+1,0,1:,5].squeeze()

# degree 200
y1_loading_degree200 = kernels[degrees[2]+1,0,1:,0].squeeze()
y2_loading_degree200 = kernels[degrees[2]+1,0,1:,1].squeeze()
y3_loading_degree200 = kernels[degrees[2]+1,0,1:,2].squeeze()
y4_loading_degree200 = kernels[degrees[2]+1,0,1:,3].squeeze()
y5_loading_degree200 = kernels[degrees[2]+1,0,1:,4].squeeze()
y6_loading_degree200 = kernels[degrees[2]+1,0,1:,5].squeeze()
#}}}

# Validate loading solutions against the analytic solutions
if validation:
    param['source'] = md.love.forcing_type

    # Extract analytic solutions
    for jj in range(0, len(degrees_analytic)):
        param['degree'] = degrees_analytic[jj]
        y_temp = yi_analytic_homogenous(param)
        if 'y_ana' not in locals(): # NOTE: Probably not a very Pythonic thing to do, but translated from MATLAB
            y_ana = np.zeros((len(degrees_analytic), np.shape(y_temp)[0], np.shape(y_temp)[1]))
        y_ana[jj,:,:] = y_temp

    y1_ana = np.squeeze(y_ana[:,:,0])
    y2_ana = np.squeeze(y_ana[:,:,1])
    y3_ana = np.squeeze(y_ana[:,:,2])
    y4_ana = np.squeeze(y_ana[:,:,3])
    y5_ana = np.squeeze(y_ana[:,:,4])
    y6_ana = np.squeeze(y_ana[:,:,5])

    depth = (np.max(param['radius']) - param['radius']) / 1000 # km

    N = 6 * (md.materials.numlayers + 1)
    kernels = np.reshape(md.results.LoveSolution.LoveKernels, [N, md.love.sh_nmax + 1, md.love.nfreq])
    kernels = kernels[[*range(3, N), *range(0, 2 + 1)], :, :]

    kernels = np.reshape(kernels, [6, (md.materials.numlayers + 1), (md.love.sh_nmax + 1), md.love.nfreq], order='F')

    kernels[:,-1,:,0] = kernels[[0, 4, 1, 3, 2, 5], -1, :, 0]

    a = md.materials.radius[-1]
    g0 = md.love.g0
    mu0 = md.love.mu0

    y1 = np.squeeze(kernels[0,:,:,0]).T
    y2 = np.squeeze(kernels[1,:,:,0]).T
    y3 = np.squeeze(kernels[2,:,:,0]).T
    y4 = np.squeeze(kernels[3,:,:,0]).T
    y5 = np.squeeze(kernels[4,:,:,0]).T
    y6 = np.squeeze(kernels[5,:,:,0]).T

    y1 = y1 * a
    y2 = y2 * mu0
    y3 = y3 * a
    y4 = y4 * mu0
    y5 = y5 * g0 * a
    y6 = y6 * g0

    # TODO: Generate plots similar to those created in MATLAB test

field_names = [
    'y1_tidal_degree002', 'y2_tidal_degree002', 'y3_tidal_degree002', 'y4_tidal_degree002', 'y5_tidal_degree002', 'y6_tidal_degree002',
    'y1_tidal_degree020', 'y2_tidal_degree020', 'y3_tidal_degree020', 'y4_tidal_degree020', 'y5_tidal_degree020', 'y6_tidal_degree020',
    'y1_tidal_degree200', 'y2_tidal_degree200', 'y3_tidal_degree200', 'y4_tidal_degree200', 'y5_tidal_degree200', 'y6_tidal_degree200',
    'y1_loading_degree002', 'y2_loading_degree002', 'y3_loading_degree002', 'y4_loading_degree002', 'y5_loading_degree002', 'y6_loading_degree002',
    'y1_loading_degree020', 'y2_loading_degree020', 'y3_loading_degree020', 'y4_loading_degree020', 'y5_loading_degree020', 'y6_loading_degree020',
    'y1_loading_degree200', 'y2_loading_degree200', 'y3_loading_degree200', 'y4_loading_degree200', 'y5_loading_degree200', 'y6_loading_degree200'
    ]
field_tolerances = [
    3e-7, 3e-7, 3e-7, 1e-7, 6e-8, 9e-7,
    2e-7, 7e-8, 3e-7, 9e-8, 9e-10, 8e-10,
    2e-8, 4e-8, 4e-7, 3e-8, 2e-10, 2e-7,
    4e-6, 1e-6, 4e-6, 3e-6, 8e-7, 2e-7,
    2e-6, 1e-7, 5e-6, 3e-7, 2e-7, 2e-7,
    2e-6, 9e-10, 5e-5, 3e-8, 5e-7, 4e-7
    ]
field_values = [
    y1_tidal_degree002, y2_tidal_degree002, y3_tidal_degree002, y4_tidal_degree002, y5_tidal_degree002, y6_tidal_degree002,
    y1_tidal_degree020, y2_tidal_degree020, y3_tidal_degree020, y4_tidal_degree020, y5_tidal_degree020, y6_tidal_degree020,
    y1_tidal_degree200, y2_tidal_degree200, y3_tidal_degree200, y4_tidal_degree200, y5_tidal_degree200, y6_tidal_degree200,
    y1_loading_degree002, y2_loading_degree002, y3_loading_degree002, y4_loading_degree002, y5_loading_degree002, y6_loading_degree002,
    y1_loading_degree020, y2_loading_degree020, y3_loading_degree020, y4_loading_degree020, y5_loading_degree020, y6_loading_degree020,
    y1_loading_degree200, y2_loading_degree200, y3_loading_degree200, y4_loading_degree200, y5_loading_degree200, y6_loading_degree200
    ]
