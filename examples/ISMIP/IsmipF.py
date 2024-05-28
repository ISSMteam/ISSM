import numpy as np
from plotmodel import plotmodel
from SetIceSheetBC import SetIceSheetBC
#Parameterization for ISMIP F experiment

#Set the Simulation generic name #md.miscellaneous
#->

#Geometry
print('   Constructing Geometry')

#Define the geometry of the simulation #md.geometry
#surface is [-x*tan(3.0*pi/180)] #md.mesh
#->

#base is [surface-1000+100*exp(-((x-L/2).^2+(y-L/2).^2)/(10000.^2))]
#L is the size of the side of the square #max(md.mesh.x)-min(md.mesh.x)
#->

#->

#thickness is the difference between surface and base #md.geometry
#->

#plot the geometry to check it out
#->


print('   Defining friction parameters')

#These parameters will not be used but need to be fixed #md.friction
#one friciton coefficient per node (md.mesh.numberofvertices,1)
#conversion from year to seconds with #md.constants.yts
#->

#one friction exponent (p,q) per element
#->

#->


print('   Construct ice rheological properties')

#The rheology parameters sit in the material section #md.materials
#B has one value per vertex
#->

#n has one value per element
#->


print('   Set boundary conditions')

#Set the default boundary conditions for an ice-sheet
# #help SetIceSheetBC
#->

print('   Initializing velocity and pressure')

#initialize the velocity and pressurefields of #md.initialization
#->

#->

#->

#->
