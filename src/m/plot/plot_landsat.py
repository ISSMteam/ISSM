#!/usr/bin/env python3
import socket
import copy
import numpy as np
import matplotlib.pyplot as plt
from cmaptools import getcolormap
from processmesh import processmesh
from processdata import processdata
from InterpFromMeshToGrid import InterpFromMeshToGrid
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from applyoptions import applyoptions

def plot_landsat(md,data,options,fig,axgrid,gridindex):
    """
    Explain
    -------
     This funtion loads Landsat Image Mosaic Antarctica (LIMA) for background image.

    Usage
    -----
    plot_landsat(md,data,options,plotlines,plotcols,i)
    """

    #process mesh and data
    x2d, y2d, z2d, elements2d, is2d, isplanet=processmesh(md,[],options)
    data, datatype=processdata(md,data,options)

    ismask = options.exist('mask')
    if ismask:
        mask = options.getfieldvalue('mask')
        options2 = copy.deepcopy(options)
        options2.removefield('caxis',False)
        options2.removefield('log',False)
        mask, datatype=processdata(md,mask,options2)
        data[~mask] = np.nan

    #check is2d
    if not is2d:
       raise Exception('buildgridded error message: gridded not supported for 3d meshes, project on a layer')

    #Get options
    hostname = socket.gethostname().lower().replace('-','')
    transparency = options.getfieldvalue('transparency',.2)
    highres = options.getfieldvalue('highres',0)
    isunit  = options.getfieldvalue('unit',1)

    #Get xlim, and ylim
    xlim=np.array(options.getfieldvalue('xlim',[min(x2d),max(x2d)]))/isunit
    ylim=np.array(options.getfieldvalue('ylim',[min(y2d),max(y2d)]))/isunit

    pwr = md.radaroverlay.pwr
    xm  = md.radaroverlay.x
    ym  = md.radaroverlay.y
    if md.mesh.epsg == 3031 & np.size(pwr)==0 | np.size(xm)==0 | np.size(ym) == 0:
        #Antarctica region
        if highres: 
            print('   LIMA with geotiff') # {{{

            # find merged mosaic landsat image
            limapath = {'simba00':'/drive/project_inwoo/issm/Data/LIMA/AntarcticaLandsat.tif'};
            limapath = limapath[hostname]
            print('   LIMA path is %s'%(limapath))

            # read image
            #im = imread(limapath);

            ## Region of LIMA data set
            #info = gdalinfo(limapath); # get geotiff info
            #xm = info.xmin + info.dx*np.arange(info.nx)
            #ym = info.ymax - info.dy*np.arange(info.ny)

            ## find region of model at LIMA
            #offset = 1e+4;
            #posx = np.where((xm > xlim[0]-offset)&(xm < xlim[1]+offset))[0]
            #posy = np.where((ym > ylim[0]-offset)&(ym < ylim[1]+offset))[0]
            # }}}
        else:
            print('   LIMA with reduced tiff') # {{{
            #Find merged mosaic landsat image
            limapath = {'inwoob85md3h':'/drive/project_inwoo/issm/Data/LIMA/tiff_90pct/00000-20080319-092059124.tif',
                        'simba00':'/home/inwoo/data/LIMA/tiff_90pct/00000-20080319-092059124.tif'}
            if not hostname in limapath.keys():
                raise Exception('Error: Landsat image at Antarctic region is downloaded at https://lima.usgs.gov/fullcontinent.php. Download geotiff image using "wget -c https://lima.usgs.gov/tiff_90pct.zip -O tiff_90pct.zip"');

            limapath = limapath[hostname]
            print('   LIMA path is %s'%(limapath))

            # read image
            #im = imread(limapath)

            ## Region of LIMA data set
            #info = gdalinfo(limapath) # get geotiff info
            #xm = info.xmin + info.dx*np.arange(info.nx)
            #ym = info.ymax - info.dy*np.arange(info.ny)

            ## find region of model at LIMA
            #offset = 1e+4
            #posx = np.where((xm > xlim[0]-offset)&(xm < xlim[1]+offset))[0]
            #posy = np.where((ym > ylim[0]-offset)&(ym < ylim[1]+offset))[0]
            # }}}

        # Update region of radaroverlay
        md.radaroverlay.x = xm[posx]
        md.radaroverlay.y = ym[posy]
        md.radaroverlay.pwr = im[posy, posx,:]
    elif md.mesh.epsg == 3431 & np.size(pwr)==0 | np.size(xm)==0 | np.size(ym) == 0:
        #Greenland region
        raise Exception('Greenland region is not yet available.')

    #Check dataset.
    if (np.size(pwr)>0) & (np.size(xm)>0) & (np.size(ym)>0):
        #Existing radaroverlay
        if np.ndim(pwr) != 3:
            raise Exception('Error: Check np.ndim(md.radaroverlay.pwr) should be equal to 3.')

        if np.any(np.diff(xm) < 0):
            print('WARNING: md.radaroverlay.x should be increasing order.')
            xm = np.flip(xm)
            pwr= np.flip(pwr,axis=0)
        if np.any(np.diff(md.radaroverlay.y) < 0):
            print('WARNING: md.radaroverlay.y should be increasing order.')
            ym = np.flip(ym)
            pwr= np.flip(pwr,axis=1)

        #Check image size
        #shape of image should be (nx, ny, band)
        nx = len(xm)
        ny = len(ym)
        if (np.shape(pwr)[0]==nx) & (np.shape(pwr)[1]==ny):
            pwr = np.transpose(pwr,[1,0,2])

        #Close-up to specific xlim
        posx = np.where((xlim[0]<=xm)&(xm<=xlim[1]))[0]
        posy = np.where((ylim[0]<=ym)&(ym<=ylim[1]))[0]
        xm = xm[posx]
        ym = ym[posy]
        pwr = pwr[posy[0]:posy[-1]+1,posx[0]:posx[-1]+1,:]
    else:
        raise Exception('Error: data array in md.radaroverlay is not implemented yet.')

    #Prepare grid
    if np.ndim(xm) == 1:
        data_grid=InterpFromMeshToGrid(elements2d+1,x2d/isunit,y2d/isunit,data,xm,ym,np.nan)
    else:
        data_grid=InterpFromMeshToMesh2d(elements2d+1,x2d,y2d,data,np.ravel(xm),np.ravel(ym),'default',np.nan)
        data_grid=np.reshape(data_grid,np.shape(xm))

    data_nan=np.isnan(data_grid)
    if options.exist('caxis'):
       caxis_opt=options.getfieldvalue('caxis')
       data_grid[np.where(data_grid<caxis_opt[0])]=caxis_opt[0]
       data_grid[np.where(data_grid>caxis_opt[1])]=caxis_opt[1]
       data_min=caxis_opt[0];
       data_max=caxis_opt[1];
    else:
       data_min=np.nanmin(data_grid)
       data_max=np.nanmax(data_grid)
       caxis_opt=[data_min, data_max] # Back-up caxis for 'applyoptions'

    options.addfielddefault('colormap',plt.cm.viridis)
    cmap = getcolormap(copy.deepcopy(options))
    #TODO: Matlab version
    #image_rgb = ind2rgb(uint16((data_grid - data_min)*(length(colorm)/(data_max-data_min))),colorm);
    #Set log scale
    if options.exist('log'):
        #NOTE: Tricy part for log scale dataset. "log" scale option does not rely on "processdata.py" function.
        data_grid=np.log(data_grid)/np.log(options.getfieldvalue('log'))
        data_min =np.log(data_min)/np.log(options.getfieldvalue('log'))
        data_max =np.log(data_max)/np.log(options.getfieldvalue('log'))

    #NOTE: Python version for ind2rgb
    image_rgb = cmap((data_grid-data_min)/(data_max-data_min))

    alpha=np.ones(np.shape(data_grid))
    alpha[np.where(~data_nan)]=transparency
    alpha=np.repeat(alpha[:,:,np.newaxis],axis=2,repeats=3)

    final=alpha*(pwr/255)+(1-alpha)*image_rgb[:,:,:3]

    #Select plot area 
    ax = axgrid[gridindex]

    xmin = min(xm)/isunit
    xmax = max(xm)/isunit
    ymin = min(ym)/isunit
    ymax = max(ym)/isunit

    #Draw RGB image
    h=ax.imshow(final, extent=[xmin, xmax, ymin, ymax], origin='lower')

    #last step: mesh gridded?
    if options.exist('edgecolor'):
        ax.triplot(x,y,triangles=elements,
                   color=options.getfieldvalue('edgecolor'),
                   linewdith=options.getfieldvalue('linewidth',1),
                   )

    #Apply options
    if ~np.isnan(data_min):
        options.changefieldvalue('caxis',caxis_opt) # force caxis so that the colorbar is ready
    options.addfielddefault('axis','xy equal off') # default axis
    applyoptions(md,data,options,fig,axgrid,gridindex)
