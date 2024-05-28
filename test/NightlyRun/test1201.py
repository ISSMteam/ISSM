#Test Name: EISMINTMassConservation
import numpy as np
from model import *
from socket import gethostname
from bamg import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


"""
This test is a test from the EISMINT for Ice shelves Vincent Rommelaere 1996.
"""

printingflag = False
results = []

for stabilization in range(1, 4):
    #The goal is to test the masstransport model
    md = bamg(model(), 'domain', '../Exp/SquareEISMINT.exp', 'hmax', 3000.)
    md = setmask(md, 'all', '')
    md = parameterize(md, '../Par/SquareEISMINT.py')
    md.smb.mass_balance[:] = 0.
    md = setflowequation(md, 'SSA', 'all')
    md.cluster = generic('name', gethostname(), 'np', 8)

    print("      initial velocity")
    md.initialization.vx = np.zeros((md.mesh.numberofvertices))
    md.initialization.vy = -400. * np.ones((md.mesh.numberofvertices))

#Stabilization
    if stabilization == 2:
        md.masstransport.stabilization = 0
    else:
        md.masstransport.stabilization = stabilization

#spc thickness
    pos = np.where(md.mesh.y > 199999.9)[0]
    times = np.arange(0, 501)
    md.masstransport.spcthickness = np.nan * np.ones((md.mesh.numberofvertices + 1, np.size(times)))
    md.masstransport.spcthickness[-1, :] = times
    md.masstransport.spcthickness[pos, :] = 500. + 100. * np.sin(2. * np.pi * times / 200.)
    if stabilization == 3:
        pos = np.nonzero(np.isnan(md.masstransport.spcthickness))
        md.masstransport.spcthickness[pos] = 500.  #No NaN for DG

#solve
    md.transient.isstressbalance = False
    md.settings.output_frequency = 500  #keep only last step
    md.verbose = verbose()
    md = solve(md, 'Transient')
    results.append(md.results.TransientSolution[-1].Thickness)

#plot results
#[elements, x, y, z, s, h1] = SectionValues(md, results[0], '../Exp/CrossLineEISMINT.exp', 100.)
#[elements, x, y, z, s, h2] = SectionValues(md, results[1], '../Exp/CrossLineEISMINT.exp', 100.)
#[elements, x, y, z, s, h3] = SectionValues(md, results[2], '../Exp/CrossLineEISMINT.exp', 100.)
#[elements, x, y, z, s, hth] = SectionValues(md, 500 + 100 * sin(2 * pi / 200 * (500 - md.mesh.y / 400)), '../Exp/CrossLineEISMINT.exp', 100.)
#plot(s, h1, 'r', s, h2, 'b', s, h3, 'g', s, hth, 'k')
#legend('Art. diff.', 'No Art. diff.', 'D.G.', 'Theoretical')
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       export_fig([issmdir() '/website/doc_pdf/validation/Images/EISMINT/IceShelf/eismintmasscthickness.pdf'])

#Fields and tolerances to track changes
field_names = ['ThicknessArtDiff', 'ThicknessNoArtDiff', 'ThicknessDG']
field_tolerances = [1e-13, 1e-13, 1e-13]
field_values = [results[0], results[1], results[2]]
