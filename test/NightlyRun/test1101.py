#Test Name: ISMIPAHO
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
    md = parameterize(md, '../Par/ISMIPA.py')
    md.extrude(9, 1.)

    md = setflowequation(md, 'HO', 'all')

#Create dirichlet on the bed only
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

    pos = np.where(md.mesh.vertexonbase)
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.

#Create MPCs to have periodic boundary conditions this is done on matlab indexing
    posx = np.where(md.mesh.x == 0.)[0]
    posx2 = np.where(md.mesh.x == np.max(md.mesh.x))[0]

    posy = np.where(np.logical_and.reduce((md.mesh.y == 0., md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]  #Don't take the same nodes two times
    posy2 = np.where(np.logical_and.reduce((md.mesh.y == np.max(md.mesh.y), md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]

    md.stressbalance.vertex_pairing = np.vstack((np.vstack((posx + 1, posx2 + 1)).T, np.vstack((posy + 1, posy2 + 1)).T))

#Compute the stressbalance
    md.cluster = generic('name', gethostname(), 'np', 8)
    md = solve(md, 'Stressbalance')

#Plot the results and save them
    vx = md.results.StressbalanceSolution.Vx
    vy = md.results.StressbalanceSolution.Vy
    vz = md.results.StressbalanceSolution.Vz
    results.append(md.results.StressbalanceSolution)
    minvx.append(np.min(vx[-md.mesh.numberofvertices2d:]))
    maxvx.append(np.max(vx[-md.mesh.numberofvertices2d:]))

#Now plot vx, vy, vz and vx on a cross section
#   plotmodel(md, 'data', vx, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km')
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipaHOvx' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipaHOvx%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
#   plotmodel(md, 'data', vy, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km')
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipaHOvy' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipaHOvy%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
#   plotmodel(md, 'data', vz, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km')
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipaHOvz' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipaHOvz%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

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
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipaHOvxsec' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipaHOvxsec%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

#Now plot the min and max values of vx for each size of the square
#plot([5 10 20 40 80 160], minvx)ylim([0 18])xlim([0 160])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipaHOminvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipaHOminvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')
#plot([5 10 20 40 80 160], maxvx)ylim([0 120])xlim([0 160])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipaHOmaxvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipaHOmaxvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestA')

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-08, 1e-08, 1e-08]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
