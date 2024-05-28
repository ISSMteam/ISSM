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

# For volumetric potential
md = model()
md.groundingline.migration = 'None'

md.materials = materials('litho')
cst = 365.25 * 24 * 3600 * 1000

md.materials.numlayers = 40
md.love.forcing_type = 9

md.materials.density = np.zeros((md.materials.numlayers, 1)) + 5511
md.materials.lame_mu = np.zeros((md.materials.numlayers, 1)) + 0.75e11
md.materials.lame_lambda = np.zeros((md.materials.numlayers, 1)) + 5e17
md.materials.issolid = np.ones((md.materials.numlayers, 1))
md.materials.rheologymodel = np.zeros((md.materials.numlayers, 1))

# The following isn't used here but needs to hhave arrays of consistent size with the rest of the materials
md.materials.viscosity = np.zeros((md.materials.numlayers, 1)) + 1e21
md.materials.burgers_mu = md.materials.lame_mu / 3
md.materials.burgers_viscosity = md.materials.viscosity / 10
md.materials.ebm_alpha = np.ones((md.materials.numlayers, 1)) * 0.9
md.materials.ebm_delta = np.ones((md.materials.numlayers, 1)) * 0.2
md.materials.ebm_taul = np.ones((md.materials.numlayers, 1)) * 54 * 60 # 54 min
md.materials.ebm_tauh = np.ones((md.materials.numlayers, 1)) * 18.6 * cst / 1e3 # 18.6 yr

md.materials.radius = np.linspace(10e3, 6371e3, md.materials.numlayers + 1).reshape(-1, 1)
md.love.g0 = 9.8134357285509388 # directly grabbed from fourierlovesolver for this particular case

md.love.allow_layer_deletion = 1
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

# For surface load
md.love.forcing_type = 11
md = solve(md,'lv')
kernels = np.reshape(md.results.LoveSolution.LoveKernels, (md.love.sh_nmax + 1, md.love.nfreq, md.materials.numlayers + 1, 6), order='F')

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
    2e-8, 4e-8, 4e-7, 3e-8, 2e-10, 1e-10,
    4e-6, 1e-6, 4e-6, 3e-6, 8e-7, 2e-6,
    2e-6, 1e-7, 5e-6, 3e-7, 2e-7, 2e-7,
    2e-6, 9e-10, 5e-5, 3e-8, 5e-7, 2e-9
    ]
field_values = [
    y1_tidal_degree002, y2_tidal_degree002, y3_tidal_degree002, y4_tidal_degree002, y5_tidal_degree002, y6_tidal_degree002,
    y1_tidal_degree020, y2_tidal_degree020, y3_tidal_degree020, y4_tidal_degree020, y5_tidal_degree020, y6_tidal_degree020,
    y1_tidal_degree200, y2_tidal_degree200, y3_tidal_degree200, y4_tidal_degree200, y5_tidal_degree200, y6_tidal_degree200,
    y1_loading_degree002, y2_loading_degree002, y3_loading_degree002, y4_loading_degree002, y5_loading_degree002, y6_loading_degree002,
    y1_loading_degree020, y2_loading_degree020, y3_loading_degree020, y4_loading_degree020, y5_loading_degree020, y6_loading_degree020,
    y1_loading_degree200, y2_loading_degree200, y3_loading_degree200, y4_loading_degree200, y5_loading_degree200, y6_loading_degree200
    ]
