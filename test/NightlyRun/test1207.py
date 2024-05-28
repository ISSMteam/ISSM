#Test Name: EISMINTRoundIceSheetStaticFS
import numpy as np
from model import *
from socket import gethostname
from roundmesh import *
from setmask import *
from parameterize import *
from setflowequation import *
from solve import *


"""
The aim of this program is to compare a model with an analytical solution given in SSA EISMINT : Lessons in Ice-Sheet Modeling.
"""

printingflag = False

numlayers = 10
resolution = 30000.

#To begin with the numerical model
md = model()
md = roundmesh(md, 750000., resolution)
md = setmask(md, '', '')  #We can not test iceshelves nor ice rises with this analytical solution
md = parameterize(md, '../Par/RoundSheetStaticEISMINT.py')

#Calculation of the analytical 2d velocity field
constant = 0.3
vx_obs = constant / 2. * md.mesh.x * (md.geometry.thickness)**- 1
vy_obs = constant / 2. * md.mesh.y * (md.geometry.thickness)**- 1
vel_obs = np.sqrt((md.inversion.vx_obs)**2 + (md.inversion.vy_obs)**2)

#We extrude the model to have a 3d model
md.extrude(numlayers, 1.)
md = setflowequation(md, 'FS', 'all')

#Spc the nodes on the bed
pos = np.where(md.mesh.vertexonbase)
md.stressbalance.spcvx[pos] = 0.
md.stressbalance.spcvy[pos] = 0.
md.stressbalance.spcvz[pos] = 0.

#Now we can solve the problem
md.cluster = generic('name', gethostname(), 'np', 8)
md = solve(md, 'Stressbalance')

#Calculate the depth averaged velocity field (2d):
vx = md.results.StressbalanceSolution.Vx
vy = md.results.StressbalanceSolution.Vy
vel = np.zeros((md.mesh.numberofvertices2d))

for i in range(0, md.mesh.numberofvertices2d):
    node_vel = 0.
    for j in range(0, md.mesh.numberoflayers - 1):
        node_vel = node_vel + (1. / (2. * (md.mesh.numberoflayers - 1)) * (np.sqrt(vx[i + (j + 1) * md.mesh.numberofvertices2d, 0]**2 + vy[i + (j + 1) * md.mesh.numberofvertices2d, 0]**2) + np.sqrt(vx[i + j * md.mesh.numberofvertices2d, 0]**2 + vy[i + j * md.mesh.numberofvertices2d, 0]**2)))
    vel[i] = node_vel

#Plot of the velocity from the exact and calculated solutions
#figure(1)
#subplot(2, 2, 1)
#p = patch('Faces', md.mesh.elements2d, 'Vertices', [md.mesh.x2d md.mesh.y2d], 'FaceVertexCData', ...
#vel, 'FaceColor', 'interp', 'EdgeColor', 'none')
#title('Modelled velocity', 'FontSize', 14, 'FontWeight', 'bold')
#colorbar
#caxis([0 200])

#subplot(2, 2, 2)
#p = patch('Faces', md.mesh.elements2d, 'Vertices', [md.mesh.x2d md.mesh.y2d], 'FaceVertexCData', ...
#vel_obs, 'FaceColor', 'interp', 'EdgeColor', 'none')
#title('Analytical velocity', 'FontSize', 14, 'FontWeight', 'bold')
#colorbar
#caxis([0 200])

#subplot(2, 2, 3)
#hold on
#plot(sqrt((md.mesh.x(1:md.mesh.numberofvertices2d)).^2 + (md.mesh.y(1:md.mesh.numberofvertices2d)).^2), vel, 'r.')
#plot(sqrt((md.mesh.x2d).^2 + (md.mesh.y2d).^2), vel_obs, 'b.')
#title('Analytical vs calculated velocity', 'FontSize', 14, 'FontWeight', 'bold')
#xlabel('distance to the center of the icesheet [m]', 'FontSize', 14, 'FontWeight', 'bold')
#ylabel('velocity [m / yr]', 'FontSize', 14, 'FontWeight', 'bold')
#legend('calculated velocity', 'exact velocity')
#axis([0 750000 0 200])
#hold off

#subplot(2, 2, 4)
#p = patch('Faces', md.mesh.elements2d, 'Vertices', [md.mesh.x2d md.mesh.y2d], 'FaceVertexCData', ...
#abs(vel - vel_obs). / vel_obs * 100, 'FaceColor', 'interp', 'EdgeColor', 'none')
#title('Relative misfit [%]', 'FontSize', 14, 'FontWeight', 'bold')
#colorbar
#caxis([0 100])

if printingflag:
    pass
#       set(gcf, 'Color', 'w')
#       printmodel('FSstatic', 'png', 'margin', 'on', 'marginsize', 25, 'frame', 'off', 'resolution', 0.7, 'hardcopy', 'off')
#       system(['mv FSstatic.png ' ISSM_DIR '/website/doc_pdf/validation/Images/EISMINT/IceSheet'])

#Fields and tolerances to track changes
field_names = ['Vx', 'Vy', 'Vel']
field_tolerances = [1e-12, 1e-12, 1e-12]
field_values = [vx, vy, vel]
