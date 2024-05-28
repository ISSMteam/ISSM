#!/usr/bin/env python
from model import model
import numpy as np
from GetAreas import GetAreas
from mesh3dprisms import mesh3dprisms
from mesh2d import mesh2d

def VolumeAboveFloatation(md, step=None, flags=None):
   """VolumeAboveFlotation - returns volume above floatation

   Usage:
      V = VolumeAboveFloatation(md)          % uses model fiels alone
      V = VolumeAboveFloatation(md,10)       % Will look at step 10 of transient solution
      V = VolumeAboveFloatation(md,10,flags) % Will look at step 10 of transient solution, only flaged elements
   """
   isverb = 0 # verbosity.

   #Special case if 3d
   if isinstance(md.mesh, mesh3dprisms):
      index = md.mesh.elements2d-1;
      x = md.mesh.x2d;
      y = md.mesh.y2d;
   elif isinstance(md.mesh, mesh2d):
      index = md.mesh.elements-1;
      x = md.mesh.x;
      y = md.mesh.y;
   else:
      raise Exception('not supported yet for {}.'%(type(md.mesh)));

   #1. get some parameters
   rho_ice   = md.materials.rho_ice
   rho_water = md.materials.rho_water

   #2. compute averages
   if (not step) and (not flags):
      base           = np.mean(md.geometry.base[index],axis=1);
      surface        = np.mean(md.geometry.surface[index],axis=1);
      bathymetry     = np.mean(md.geometry.bed[index],axis=1);
      ice_levelset   = md.mask.ice_levelset;
      ocean_levelset = md.mask.ocean_levelset;
   else:
      if 'MaskIceLevelset' in md.results.TransientSolution[step].keys():
      #if isprop(md.results.TransientSolution(step),'MaskIceLevelset')
         ice_levelset   = md.results.TransientSolution[step].MaskIceLevelset;
      else:
         ice_levelset   = md.mask.ice_levelset;
      ocean_levelset = md.results.TransientSolution[step].MaskOceanLevelset;
      base           = np.mean(md.results.TransientSolution[step].Base[index],axis=1);
      surface        = np.mean(md.results.TransientSolution[step].Surface[index],axis=1);
      if 'Bed' in md.results.TransientSolution[step].keys(): #,'Bed')
         bathymetry  = np.mean(md.results.TransientSolution[step].Bed[index],axis=1);
      else:
         bathymetry  = np.mean(md.geometry.bed[index],axis=1);

   #3. get areas of all triangles
   areas = GetAreas(index+1,x,y);

   #4. Compute volume above floatation
   if isverb:
      print(np.shape(areas))
      print(np.shape(surface))
      print(np.shape(base))
      print(np.shape(bathymetry))

   V = areas*(surface-base+np.minimum(rho_water/rho_ice*bathymetry,0.))
   if isverb:
      print(np.shape(V))

   #5. take out the ones that are outside of levelset or floating
   pos = np.where((np.min(ice_levelset[index],axis=1)>0) | (np.min(ocean_levelset[index],axis=1)<0))
   V[pos] = 0;

   #In case we are only looking at one portion of the domain...
   if flags:
      V[~flags] = 0;

   #sum individual contributions
   V = np.sum(V);

   if isverb:
      print('   potential volume is: %e m^3'%(V))

   return V
