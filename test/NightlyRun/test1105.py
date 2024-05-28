#Test Name: ISMIPCHO
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

for L in L_list:  #in m (3 times the desired length for BC problems)
    nx = 30  #number of nodes in x direction
    ny = 30
    md = model()
    md = squaremesh(md, L, L, nx, ny)
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPC.py')
    md.extrude(10, 1.)

    md = setflowequation(md, 'HO', 'all')

#Create MPCs to have periodic boundary conditions
    md.stressbalance.spcvx = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvy = np.nan * np.ones((md.mesh.numberofvertices))
    md.stressbalance.spcvz = np.nan * np.ones((md.mesh.numberofvertices))

    posx = np.where(np.logical_and.reduce((md.mesh.x == 0., md.mesh.y != 0., md.mesh.y != L)))[0]
    posx2 = np.where(np.logical_and.reduce((md.mesh.x == L, md.mesh.y != 0., md.mesh.y != L)))[0]

    posy = np.where(np.logical_and.reduce((md.mesh.y == 0., md.mesh.x != 0., md.mesh.x != L)))[0]  #Don't take the same nodes two times
    posy2 = np.where(np.logical_and.reduce((md.mesh.y == L, md.mesh.x != 0., md.mesh.x != L)))[0]

    md.stressbalance.vertex_pairing = np.vstack((np.vstack((posx + 1, posx2 + 1)).T, np.vstack((posy + 1, posy2 + 1)).T))

#Add spc on the corners
    pos = np.where(np.logical_and.reduce((np.logical_or(md.mesh.x == 0., md.mesh.x == L), np.logical_or(md.mesh.y == 0., md.mesh.y == L), md.mesh.vertexonbase)))
    md.stressbalance.spcvx[pos] = 0.
    md.stressbalance.spcvy[pos] = 0.
    if (L == 5000.):
        md.stressbalance.spcvx[pos] = 15.66
        md.stressbalance.spcvy[pos] = -0.1967
    elif (L == 10000.):
        md.stressbalance.spcvx[pos] = 16.04
        md.stressbalance.spcvy[pos] = -0.1977
    elif (L == 20000.):
        md.stressbalance.spcvx[pos] = 16.53
        md.stressbalance.spcvy[pos] = -1.27
    elif (L == 40000.):
        md.stressbalance.spcvx[pos] = 17.23
        md.stressbalance.spcvy[pos] = -3.17
    elif (L == 80000.):
        md.stressbalance.spcvx[pos] = 16.68
        md.stressbalance.spcvy[pos] = -2.69
    elif (L == 160000.):
        md.stressbalance.spcvx[pos] = 16.03
        md.stressbalance.spcvy[pos] = -1.27

#Spc the bed at zero for vz
    pos = np.where(md.mesh.vertexonbase)
    md.stressbalance.spcvz[pos] = 0.

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
#   plotmodel(md, 'data', vx, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 2)
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipcHOvx' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipcHOvx%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')
#   plotmodel(md, 'data', vy, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 3)
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipcHOvy' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipcHOvy%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')
#   plotmodel(md, 'data', vz, 'layer  #all', md.mesh.numberoflayers, 'xlim', [0 L / 1.0e3], 'ylim', [0 L / 1.0e3], 'unit', 'km', 'figure', 4)
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipcHOvz' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipcHOvz%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')

    if (L == 5000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP5000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [0 20], 'xlim', [0 5000], 'title', '', 'xlabel', '', 'figure', 5)
    elif (L == 10000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP10000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [13 18], 'xlim', [0 10000], 'title', '', 'xlabel', '')
    elif (L == 20000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP20000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [14 22], 'xlim', [0 20000], 'title', '', 'xlabel', '')
    elif (L == 40000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP40000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [10 40], 'xlim', [0 40000], 'title', '', 'xlabel', '')
    elif (L == 80000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP80000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [0 80], 'xlim', [0 80000], 'title', '', 'xlabel', '')
    elif (L == 160000.):
        pass
#           plotmodel(md, 'data', vx, 'sectionvalue', '../Exp/ISMIP160000.exp', 'layer', md.mesh.numberoflayers, ...
#                   'resolution', [10 10], 'ylim', [0 200], 'xlim', [0 160000], 'title', '', 'xlabel', '')
    if printingflag:
        pass
#           set(gcf, 'Color', 'w')
#           printmodel(['ismipcHOvxsec' num2str(L)], 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#           shutil.move("ismipcHOvxsec%d.png" % L, ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')

#Now plot the min and max values of vx for each size of the square
#plot([5 10 20 40 80 160], minvx)ylim([4 18])xlim([0 160])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipcHOminvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipcHOminvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')
#plot([5 10 20 40 80 160], maxvx)ylim([0 200]) xlim([0 160])
if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('ismipcHOmaxvx', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
#       shutil.move('ismipcHOmaxvx.png', ISSM_DIR + '/website/doc_pdf/validation/Images/ISMIP/TestC')

#Fields and tolerances to track changes
field_names = ['Vx80km', 'Vy80km', 'Vz80km']
field_tolerances = [1e-09, 1e-08, 1e-08]
field_values = []
for result in results:
    field_values = field_values + [result.Vx, result.Vy, result.Vz]
