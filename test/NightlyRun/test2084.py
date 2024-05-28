#Test Name: GiaCaron
# Forward Love number solution for a viscoelastic earth, model M3-L70-V01 from
# Spada, G., Barletta, V. R., Klemann, V., Riva, R. E. M., Martinec, Z.,
# Gasperini, P., Lund, B., Wolf, D., Vermeersen, L. L. A. and King, M. A.
# (2011), A benchmark study for glacial isostatic adjustment codes. Geophysical
# Journal International, 185: 106--132. doi:10.1111/j.1365-246X.2011.04952.x

import numpy as np

from generic import generic
from materials import *
from socket import gethostname
from model import *
from solve import *

md = model()
md.cluster = generic('name', gethostname(), 'np', 8)

# Set validation=1 for comparing against the Spada benchmark
validation = 0

md.materials = materials('litho')
md.miscellaneous.name = 'FourierLoveTest'
md.groundingline.migration = 'None'

md.verbose = verbose('all')
md.verbose = verbose('1111111111111111')
cst = 365.25 * 24 * 3600 * 1000

md.materials.numlayers = 6
md.materials.radius = np.array([10, 1222.5, 3.4800e+03, 5.7010e+03, 5.9510e+03,
                                6.3010e+03, 6.3710e+03]).reshape(-1, 1) * 1e3
md.materials.density = np.array([1.0750e4, 1.0750e+04, 4.9780e+03, 3.8710e+03,
                                3.4380e+03, 3.0370e+03]).reshape(-1, 1)
md.materials.lame_mu = np.array([1e-5, 0, 2.2834e+00, 1.0549e+00, 7.0363e-01,
                                5.0605e-01]).reshape(-1, 1) * 1e11
md.materials.viscosity = np.array([0, 0, 2.0000e+00, 1.0000e+00, 1.0000e+00,
                                   1.0000e+25]).reshape(-1, 1) * 1e21
md.materials.lame_lambda = np.array(md.materials.lame_mu) * 0 + 5e17
md.materials.issolid = np.array([1, 0, 1, 1, 1, 1]).reshape(-1, 1)
md.materials.rheologymodel = np.zeros((md.materials.numlayers, 1))
md.materials.burgers_mu = md.materials.lame_mu / 3
md.materials.burgers_viscosity = md.materials.viscosity / 10
md.materials.ebm_alpha = np.ones((md.materials.numlayers, 1)) * 0.9
md.materials.ebm_delta = np.ones((md.materials.numlayers, 1)) * 0.2
md.materials.ebm_taul = np.ones((md.materials.numlayers, 1)) * 54 * 60 # 54 min
md.materials.ebm_tauh = np.ones((md.materials.numlayers, 1)) * 18.6 * cst / 1e3 # 18.6 yr
#setlitho2prem(md.materials)

md.love.allow_layer_deletion = 1
md.love.frequencies = np.vstack(([0], np.logspace(-6, 3, 1000).reshape(-1, 1) / cst))
md.love.nfreq = len(md.love.frequencies)
md.love.sh_nmin = 1
md.love.sh_nmax = 1000
md.love.underflow_tol = 1e-20
md.love.pw_threshold = 1e-3
md.love.Gravitational_Constant = 6.6732e-11
md.love.allow_layer_deletion = 1
md.love.forcing_type = 11
md.love.chandler_wobble = 0
md.love.complex_computation = 0

md.love.istemporal = 1
md.love.n_temporal_iterations = 8
#md.love.time = np.logspace(-6, 5, 2).reshape(-1, 1) * cst
md.love.time = np.vstack(([0], np.logspace(-3, 5, 24).reshape(-1, 1) * cst))

#md.love.time = np.linspace(1/12, 10, 10 * 12).reshape(-1, 1) * cst / 1e3
md.love.love_kernels = 1
if md.love.istemporal:
    md.love = md.love.build_frequencies_from_time()

md = solve(md, 'lv')

ht = md.results.LoveSolution.LoveHt.reshape(-1, 1)
lt = md.results.LoveSolution.LoveLt.reshape(-1, 1)
kt = md.results.LoveSolution.LoveKt.reshape(-1, 1)
t = md.love.time / cst * 1e3

# Fields and tolerances to track changes
# loading love numbers
field_names = ['LoveH_loading_elastic', 'LoveK_loading_elastic', 'LoveL_loading_elastic']
field_tolerances = [2.0e-8, 2.0e-8, 2.0e-8]
field_values = [
    np.array(md.results.LoveSolution.LoveHt)[:, 0],
    np.array(md.results.LoveSolution.LoveKt)[:, 0],
    np.array(md.results.LoveSolution.LoveLt)[:, 0]
]

# TODO:
# - Implement read from file and comparison
# - Implement plot
# - Implement validation of elastic loading solutions against the Spada benchmark
