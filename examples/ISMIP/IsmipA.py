import numpy as np
from plotmodel import plotmodel
from SetIceSheetBC import SetIceSheetBC
#Parameterization for ISMIP A experiment

#Set the Simulation generic name #md.miscellaneous
#->

#Geometry
print('   Constructing Geometry')

#Define the geometry of the simulation #md.geometry
#surface is [-x*tan(0.5*pi/180)] #md.mesh
#->

#base is [surface-1000+500*sin(x*2*pi/L).*sin(y*2*pi/L)]
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
#->

#one friciton exponent (p,q) per element
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
# help SetIceSheetBC
#->
