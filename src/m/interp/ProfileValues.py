#!/usr/bin/env python3
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from InterpFromMeshToMesh3d import InterpFromMeshToMesh3d
from project2d import project2d
import numpy as np
from numpy import ceil, mean, linspace, arange, transpose, reshape, ones, shape 

def ProfileValues(md,data,xprof,yprof,resolution):
    """
    PROFILEVALUES - compute the value of a field on a vertical profile

       This routine gets the value of a given field of the model on 
       a point given by its coordinates

       Usage:
          z,data=ProfileValues(md,data,xcoord,ycoord,resolution)
    """

    #Get bed and surface for each 2d point, offset to make sure that it is inside the glacier system
    offset=10e-3
    bed=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.base,1),xprof,yprof)+offset
    surface=InterpFromMeshToMesh2d(md.mesh.elements2d,md.mesh.x2d,md.mesh.y2d,project2d(md,md.geometry.surface,1),xprof,yprof)-offset

    #Some useful parameters
    layers=ceil(mean(md.geometry.thickness)/resolution)
    Z=arange(bed,surface+resolution,resolution)
    Z=transpose(Z)
    X=xprof*ones(shape(Z))
    Y=yprof*ones(shape(Z))
    data_interp=InterpFromMeshToMesh3d(md.mesh.elements,md.mesh.x,md.mesh.y,md.mesh.z,data,X,Y,Z,np.nan)

    #Get np.ndarray in list
    data_interp=data_interp[0]

    return Z, data_interp
