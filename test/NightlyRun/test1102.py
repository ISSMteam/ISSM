#Test Name: ISMIPAFS
import numpy as np
from model import *
from socket import gethostname
from squaremesh import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *

"""
This test is a test from the ISMP - HOM Intercomparison project.
Pattyn and Payne 2006
"""

printingflag = False

L_list = [80000.]
results = []
minvx = []
maxvx = []

for L in L_list:
    nx = 20  #numberof nodes in x direction
    ny = 20
    md = model()
    md = squaremesh(md, L, L, nx, ny)
    md = setmask(md, '', '')  #ice sheet test

    #  #Find elements at the corner and extract model
    #   posnodes = np.nonzero(np.logical_and.reduce(np.logical_or.reduce(md.mesh.x = 0., md.mesh.x = np.max(md.mesh.x)), np.logical_or.reduce(md.mesh.y = 0., md.mesh.y = np.max(md.mesh.y))))
    #   a = np.nonzero(ismember(md.mesh.elements, posnodes))[0]
    #   elements = np.ones((md.mesh.numberofelements), int)
    #   elements[a] = 0
    #   md.modelextract(elements)

    md = parameterize(md, '../Par/ISMIPA.py')
    md.extrude(10, 1.)
    md = setflowequation(md, 'FS', 'all')

    #Create dirichlet on the bed only
    pos = np.nonzero(md.mesh.vertexonbase)
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.
    md.stressbalance.spcvz[pos] = 0.

    #Compute the stressbalance
    md.stressbalance.abstol = np.nan
    md.stressbalance.reltol = np.nan
    md.stressbalance.restol = 1.
    md.cluster = generic('name', gethostname(), 'np', 8)
    md = solve(md, 'Stressbalance')

#Plot the results and save them
    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    pressure = md.results.StressbalanceSolution.Pressure
    results.append(md.results.StressbalanceSolution)
    minvx.append(np.min(vx[-md.mesh.numberofvertices2d:]))
    maxvx.append(np.max(vx[-md.mesh.numberofvertices2d:]))

    #Now plot vx, vy, vz and vx on a cross section
    #   plotmodel(md, 'data', vx, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 2)
    if printingflag:
        pass
    #           set(gcf, 'Color', 'w')
    #           printmodel(['ismipaFSvx' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
    #           shutil.move("ismipaFSvx%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
    #   plotmodel(md, 'data', vy, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 3)
    if printingflag:
        pass
    #           set(gcf, 'Color', 'w')
    #           printmodel(['ismipaFSvy' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
    #           shutil.move("ismipaFSvy%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
    #   plotmodel(md, 'data', vz, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 4)
    if printingflag:
        pass
    #           set(gcf, 'Color', 'w')
    #           printmodel(['ismipaFSvz' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
    #           shutil.move("ismipaFSvz%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

    if (L == 5000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP5000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [10 18], 'xlim', [0 5000], 'title', '', 'xlabel', '')
    elif (L == 10000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP10000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [10 30], 'xlim', [0 10000], 'title', '', 'xlabel', '')
    elif (L == 20000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP20000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [0 50], 'xlim', [0 20000], 'title', '', 'xlabel', '')
    elif (L == 40000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP40000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [0 80], 'xlim', [0 40000], 'title', '', 'xlabel', '')
    elif (L == 80000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP80000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [0 100], 'xlim', [0 80000], 'title', '', 'xlabel', '')
    elif (L == 160000.):
        pass
    #           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP160000.exp', 'layer', md.mesh.numberoflayers, ...
    #                   'resolution', [10 10], 'ylim', [0 120], 'xlim', [0 160000], 'title', '', 'xlabel', '')
    if printingflag:
        pass
    # set(gcf, 'Color', 'w')
    # printmodel(['ismipaFSvxsec' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
    # shutil.move("ismipaFSvxsec.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

    #Now plot the min and max values of vx for each size of the square
    #plot([5 10 20 40 80 160], minvx)ylim([0 18])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipaFSminvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipaFSminvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
#plot([5 10 20 40 80 160], maxvx)ylim([0 120])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipaFSmaxvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipaFSmaxvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-12, 1e-11, 1e-12]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
