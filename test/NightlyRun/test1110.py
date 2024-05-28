#Test Name: ISMIPF
import numpy as np
from model import *
from socket import gethostname
from bamg import *
from setmask import setmask
from parameterize import parameterize
from setflowequation import setflowequation
from solve import solve
from squaremesh import squaremesh

#This test is a test from the ISMP - HOM Intercomparison project.
#TestF
printingflag = False
results = []

for i in [1]:  #range(4):
    L = 100000.  #in m
    nx = 30  #numberof nodes in x direction
    ny = 30
    md = model()
    md = squaremesh(md, L, L, nx, ny)
    #   md = triangle(md, '../Exp/SquareISMIP.exp', 5500.)
    md = setmask(md, '', '')  #ice sheet test
    md = parameterize(md, '../Par/ISMIPF.py')
    md = md.extrude(4, 1.)
    if (i == 0 or i == 1):
        md = setflowequation(md, 'HO', 'all')
    else:
        md = setflowequation(md, 'FS', 'all')

    md.stressbalance.spcvx = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    md.stressbalance.spcvy = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    md.stressbalance.spcvz = float('NaN') * np.ones((md.mesh.numberofvertices, ))
    if (i == 0 or i == 2):
        #Create dirichlet on the bed if no slip
        pos = np.where(md.mesh.vertexonbase)
        md.stressbalance.spcvx[pos] = 0.
        md.stressbalance.spcvy[pos] = 0.
        md.stressbalance.spcvz[pos] = 0.
    else:
        posA = np.where(md.mesh.vertexonbase)
        posB = np.unique(np.hstack((np.where(md.mesh.x == 0.), np.where(md.mesh.x == np.nanmax(md.mesh.x)))))
        posC = np.unique(np.hstack((np.where(md.mesh.y == 0.), np.where(md.mesh.y == np.nanmax(md.mesh.y)))))
        pos = np.intersect1d(np.intersect1d(posA, posB), posC)
        md.stressbalance.spcvx[pos] = 100.  #because we need a dirichlet somewhere
        md.stressbalance.spcvy[pos] = 0.
        md.stressbalance.spcvz[pos] = 0.

    pos = np.where(np.logical_not(md.mesh.vertexonbase))
    md.thermal.spctemperature[pos] = 255.

    #Create MPCs to have periodic boundary conditions
    posx = np.where(md.mesh.x == 0.)[0]
    posx2 = np.where(md.mesh.x == max(md.mesh.x))[0]
    # posy = np.where(np.logical_and.reduce((md.mesh.y == 0., md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]  #Don't take the same nodes two times
    # posy2 = np.where(np.logical_and.reduce((md.mesh.y == np.max(md.mesh.y), md.mesh.x != 0., md.mesh.x != np.max(md.mesh.x))))[0]
    posy = np.where(md.mesh.y == 0)[0]
    posy2 = np.where(md.mesh.y == np.max(md.mesh.y))[0]

    md.stressbalance.vertex_pairing = np.vstack((np.vstack((posx + 1, posx2 + 1)).T, np.vstack((posy + 1, posy2 + 1)).T))
    md.masstransport.vertex_pairing = np.vstack((np.vstack((posx + 1, posx2 + 1)).T, np.vstack((posy + 1, posy2 + 1)).T))

    md.timestepping.time_step = 3.
    md.timestepping.final_time = 30.  #300.
    md.settings.output_frequency = 5  #50
    md.masstransport.stabilization = 1
    md.stressbalance.maxiter = 1

    #Compute the stressbalance
    md.cluster = generic('name', gethostname(), 'np', 8)
    md.verbose = verbose('convergence', True, 'solution', True)
    md = solve(md, 'Transient')

    #save the results
    results = np.append(results, md.results.TransientSolution[-1])

    #Now plot vx and delta surface
    if (i == 0 or i == 2):
        plotmodel(md, 'data', np.squeeze(md.results.TransientSolution[-1].Vx),
                  'layer', md.mesh.numberoflayers,
                  'sectionvalue', '../Exp/ISMIP100000.exp',
                  'title', '',
                  'xlabel', '',
                  'ylabel', 'Velocity (m / yr)',
                  'linewidth', 3,
                  'grid', 'on',
                  'unit', 'km',
                  'ylim', [91, 100])
    elif (i == 1 or i == 3):
        plotmodel(md, 'data', np.squeeze(md.results.TransientSolution[-1].Vx),
                  'layer', md.mesh.numberoflayers,
                  'sectionvalue', '../Exp/ISMIP100000.exp',
                  'title', '',
                  'xlabel', '',
                  'ylabel', 'Velocity (m / yr)',
                  'linewidth', 3,
                  'grid', 'on',
                  'unit', 'km',
                  'ylim', [185, 200])

    if printingflag:
        #set(gcf, 'Color', 'w')
        if i == 0:
            printmodel('ismipfHOvxfrozen', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfHOvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 1:
            printmodel('ismipfHOvxsliding', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfHOvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 2:
            printmodel('ismipfFSvxfrozen', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfFSvxfrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 3:
            printmodel('ismipfFSvxsliding', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfFSvxsliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])

    plotmodel(md, 'data', np.squeeze(md.results.TransientSolution[-1].Surface) - md.geometry.surface,
              'layer', md.mesh.numberoflayers,
              'sectionvalue', '../Exp/ISMIP100000.exp',
              'title', '',
              'xlabel', '',
              'ylabel', 'Surface (m)',
              'linewidth', 3,
              'grid', 'on',
              'unit', 'km',
              'ylim', [- 30, 50])
    if printingflag:
        #set(gcf, 'Color', 'w')
        if i == 0:
            printmodel('ismipfHOdeltasurfacefrozen', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfHOdeltasurfacefrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 1:
            printmodel('ismipfHOdeltasurfacesliding', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfHOdeltasurfacesliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 2:
            printmodel('ismipfFSdeltasurfacefrozen', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfFSdeltasurfacefrozen.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])
        elif i == 3:
            printmodel('ismipfFSdeltasurfacesliding', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 1.5, 'hardcopy', 'off')
            #system(['mv ismipfFSdeltasurfacesliding.png ' ISSM_DIR '/website/doc_pdf/validation/Images/ISMIP/TestF'])

#Fields and tolerances to track changes
field_names = ['VxPattynFrozen', 'VyPattynFrozen', 'VzPattynFrozen', 'SurfacePattynFrozen',
               'VxPattynSliding', 'VyPattynSliding', 'VzPattynSliding', 'SurfacePattynSliding',
               'VxStokesFrozen', 'VyStokesFrozen', 'VzStokesFrozen', 'SurfaceStokesFrozen',
               'VxStokesSliding', 'VyStokesSliding', 'VzStokesSliding', 'SurfaceStokesSliding']
field_tolerances = [1e-10, 1e-09, 1e-09, 1e-10,
                    1e-10, 1e-09, 1e-09, 1e-10,
                    1e-08, 1e-09, 1e-08, 1e-09,
                    1e-08, 2e-09, 1e-08, 1e-09]
field_values = []
for i in range(4):
    result = results[i]
    field_values += ([result.Vx, result.Vy, result.Vz, result.Surface] - md.geometry.surface)
