#Test Name: GiaIvinsBenchmarksAB2dD
from socket import gethostname

import numpy as np

from model import *
from parameterize import *
from setmask import *
from solve import *
from triangle import *


# Benchmark experiments (Figure A2a Ivins and James, 1999, Geophys. J. Int.)
md = triangle(model(), '../Exp/RoundFrontEISMINT.exp', 200000)
md = setmask(md, '', '')
md = parameterize(md, '../Par/GiaIvinsBenchmarksCD.py')

# indicate what you want to compute
md.gia.cross_section_shape = 2  # for elliptical edge

# evaluation time (termed start_time)
md.timestepping.start_time = 0.3  # for t \approx 0 kyr : to get elastic response!
md.timestepping.final_time = 2500000  # 2,500 kyr

# define loading history
md.geometry.thickness = np.array([
    np.append(md.geometry.thickness * 0.0, 0.0),
    np.append(md.geometry.thickness / 2.0, 0.1),
    np.append(md.geometry.thickness, 0.2),
    np.append(md.geometry.thickness, md.timestepping.start_time)
    ]).T

# find out elements that have zero loads throughout the loading history
pos = np.where(np.abs(md.geometry.thickness[0:-2, :].sum(axis=1)) == 0)[0]
md.mask.ice_levelset[pos] = 1 # no ice

md.cluster = generic('name', gethostname(), 'np', 3)
md.verbose = verbose('1111111')

# solve for GIA deflection
md = solve(md, 'Gia')

# Test Name: GiaIvinsBenchmarksAB2dD1
U_AB2dD1 = md.results.GiaSolution.UGia
URate_AB2dD1 = md.results.GiaSolution.UGiaRate

# Test Name: GiaIvinsBenchmarksAB2dD2
# different evaluation time # {{{
md.timestepping.start_time = 1000.3 # for t \approx 1 kyr
md.geometry.thickness[-1, -1] = md.timestepping.start_time

md = solve(md, 'Gia')

U_AB2dD2 = md.results.GiaSolution.UGia
URate_AB2dD2 = md.results.GiaSolution.UGiaRate
# }}}

# Test Name: GiaIvinsBenchmarksAB2dD3
# different evaluation time # {{{
md.timestepping.start_time = 2400000 # for t \approx \infty
md.geometry.thickness[-1, -1] = md.timestepping.start_time

md = solve(md, 'Gia')

U_AB2dD3 = md.results.GiaSolution.UGia
URate_AB2dD3 = md.results.GiaSolution.UGiaRate
# }}}

# Fields and tolerances to track changes
field_names = ['U_AB2dD1', 'URate_AB2dD1', 'U_AB2dD2', 'URate_AB2dD2', 'U_AB2dD3', 'URate_AB2dD3']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13]
field_values = [U_AB2dD1, URate_AB2dD1, U_AB2dD2, URate_AB2dD2, U_AB2dD3, URate_AB2dD3]

