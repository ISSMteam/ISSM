#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.colors
import warnings
import numpy as np
from InterpFromMeshToGrid import InterpFromMeshToGrid
from processmesh import processmesh
from processdata import processdata
from cmaptools import getcolormap
from applyoptions import applyoptions
import copy

def plot_gridded(md,data,options,fig,axgrid,gridindex):
    '''
    PLOT_OVERLAY - superimpose radar image to a given field

       Usage:
          plot_gridded(md,data,options,fig,axgrid,gridindex)

       See also: PLOTMODEL
    '''

    #process mesh and data
    x, y, z, elements, is2d, isplanet=processmesh(md,[],options)
    data, datatype=processdata(md,data,options)

    ax = axgrid[gridindex]
    #fig.delaxes(axgrid.cbar_axes[gridindex])

    islevelset = options.exist('levelset')
    if islevelset:
        levelset = options.getfieldvalue('levelset')
        options2 = copy.deepcopy(options)
        options2.removefield('caxis',False)
        options2.removefield('log',False)
        levelset, datatype=processdata(md,levelset,options2)

    #check is2d
    if not is2d:
        raise Exception('buildgridded error message: gridded not supported for 3d meshes, project on a layer')

    #Get xlim and ylim (used to extract radar image)
    xlim=options.getfieldvalue('xlim',[min(x), max(x)])
    ylim=options.getfieldvalue('ylim',[min(y), max(y)])

    isAxis = options.exist('axis')
    if isAxis:
        myaxis = options.getfieldvalue('axis');
        if isinstance(myaxis,'char'):
            xlim = [myaxis[0], myaxis[1]]
            ylim = [myaxis[2], myaxis[3]]

    postx=options.getfieldvalue('posting',np.diff(xlim)[0]/1000)
    posty=options.getfieldvalue('posting',np.diff(ylim)[0]/1000)

    #Interpolating data on grid
    nx = int(np.diff(xlim)[0]/postx)+1
    ny = int(np.diff(ylim)[0]/posty)+1
    x_m = np.linspace(xlim[0],xlim[1],nx)
    y_m = np.linspace(ylim[0],ylim[1],ny)
    #NOTE: Tricky part for elements in interpolation.
    data_grid=InterpFromMeshToGrid(elements+1,x,y,data,x_m,y_m,np.nan)
    data_grid_save = copy.deepcopy(data_grid)
    if (np.shape(data_grid)[0]<3) | (np.shape(data_grid)[1]<3):
        raise Exception('data_grid size too small in plot_gridded, check posting and units');

    #Mask values if levelset>0
    if islevelset:
        #NOTE: Tricky part for elements in interpolation.
        ls_grid=InterpFromMeshToGrid(elements+1,x,y,levelset,x_m,y_m,np.nan)
        data_grid[ls_grid>0] = np.nan

    #Process data_grid: add white in NaN and correct caxis accordingly
    data_nani, data_nanj=np.where(np.isnan(data_grid) | (data_grid==-9999))
    if options.exist('caxis'):
        caxis_opt=options.getfieldvalue('caxis');
        data_grid[np.where(data_grid<caxis_opt[0])]=caxis_opt[0]
        data_grid[np.where(data_grid>caxis_opt[1])]=caxis_opt[1]
        data_min=caxis_opt[0]
        data_max=caxis_opt[1]
    else:
        data_min=np.nanmin(data_grid[:])
        data_max=np.nanmax(data_grid[:])

    #TODO: Select plot area 
    #subplotmodel(plotlines,plotcols,i,options);

    #shading interp;
    options.addfielddefault('colormap',plt.cm.viridis)
    cmap = getcolormap(copy.deepcopy(options))
    #TODO: Matlab version
    #image_rgb = ind2rgb(uint16((data_grid - data_min)*(length(map)/(data_max-data_min))),cmap);
    #NOTE: Python version
    if isinstance(cmap,matplotlib.colors.ListedColormap):
        data_norm = (data_grid-data_min)/(data_max-data_min)
        image_rgb = cmap(data_norm)
    else:
        #TODO: Other colormaps...
        image_rgb = cmap((data_grid-data_min)/(data_max-data_min))

    #TODO: shaded...
    if options.exist('shaded'):
        warnings.warn('WARNING: shaded is not supported in Python.')

    #    if options.exist('dem'):
    #        dem, _=processdata(md,options.getfieldvalue('dem'),options)
    #        dem_grid=InterpFromMeshToGrid(elements+1,x,y,dem,x_m,y_m,np.nan);
    #    else:
    #        dem_grid=data_grid_save
    #    a    = -45;
    #    scut = 0.2;
    #    c    = 1;
    #    # computes lighting from elevation gradient
    #    fx, fy = np.gradient(dem_grid,np.gradient(x_m),np.gradient(y_m))
    #    fxy = -fx*np.sin(a*np.pi/180) - fy*np.cos(a*np.pi/180)
    #    # free some memory...
    #    del fx
    #    del fy
    #    fxy[np.isnan(fxy)] = 0

    #    # computes maximum absolute gradient (median-style), normalizes, saturates and duplicates in 3-D matrix
    #    r = np.tile(np.maximum(np.minimum(fxy/nmedian(abs(fxy),1 - scut/100),1),-1),[4,1,1])
    #    print(np.shape(r))
    #    r = np.transpose(r,[2,1,0])

    #    # applies contrast using exponent
    #    rp = (1 - abs(r))**c
    #    image_rgb = image_rgb*rp

    #    # lighter for positive gradient
    #    k = np.where(r > 0)
    #    image_rgb[k] = image_rgb[k] + (1 - rp[k])

    # set novalues / NaN to black color
    #if not np.isempty(data_nani):
    #    nancolor=options.getfieldvalue('nancolor',[1, 1, 1])
    #    image_rgb(sub2ind(size(image_rgb),repmat(data_nani,1,3),repmat(data_nanj,1,3),repmat(1:3,size(data_nani,1),1))) = repmat(nancolor,size(data_nani,1),1);

    #plot grid
    h=ax.imshow(image_rgb,extent=[xlim[0], xlim[1], ylim[0], ylim[1]],origin='lower')

    #last step: mesh gridded?
    if options.exist('edgecolor'):
        #A=elements[:,0]; B=elements[:,1]; C=elements[:,2]
        #patch('Faces',[A B C],'Vertices', [x y z],'FaceVertexCData',data_grid(1)*ones(size(x)),'FaceColor','none','EdgeColor',getfieldvalue(options,'edgecolor'));
        ax.triplot(x,y,triangles=elements,
                   color=options.getfieldvalue('edgecolor'),
                   linewdith=options.getfieldvalue('linewidth',1),
                   )

    #Apply options
    if (not np.isnan(data_min)) & (not np.isinf(data_min)):
        options.changefieldvalue('caxis',[data_min, data_max]) # force caxis so that the colorbar is ready
    options.addfielddefault('axis','xy equal'); # default axis
    applyoptions(md,data,options,fig,axgrid,gridindex)

def nmedian(x,n=0.5):
    '''
    NMEDIAN Generalized median filter
      NMEDIAN(X,N) sorts elemets of X and returns N-th value (N normalized).
      So:
         N = 0 is minimum value
         N = 0.5 is median value
         N = 1 is maximum value
    '''

    from scipy.interpolate import interp1d

    #if nargin < 2:
    #    n = 0.5;

    y = np.sort(x[:])
    xp = np.arange(0,len(y))
    f = interp1d(xp, np.sort(y), axis=0, kind='linear')
    y = f(n*(len(y)-1)+1)

    return y
