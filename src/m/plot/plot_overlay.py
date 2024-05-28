import os

import matplotlib.pyplot as plt
import matplotlib as mpl
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    print('Basemap toolkit not installed')
import numpy as np
try:
    from osgeo import gdal
except ImportError:
    print('OSGeo/GDAL for python not installed, plot_overlay is disabled')

from processdata import processdata
from processmesh import processmesh
from xy2ll import xy2ll


def plot_overlay(md, data, options, ax):
    """Function for plotting a georeferenced image. Called from plot_manager by 
    call to plotmodel.

    Usage:
        plot_overlay(md, data, options, ax)

    See also: PLOTMODEL
    """

    x, y, z, elements, is2d, isplanet = processmesh(md, [], options)
    try:
        data, datatype = processdata(md, data, options)
        imageonly = 0
    except (TypeError, ValueError):  #that should catch None and 'none' but may also catch unwanted errors
        imageonly = 1
        data = np.float('nan') * np.ones((md.mesh.numberofvertices, ))
        datatype = 1

    if not is2d:
        raise Exception('overlay plot not supported for 3D meshes, project on a 2D layer first')

    if not options.exist('overlay_image'):
        raise Exception('overlay error: provide overlay_image with path to geotiff file')
    image = options.getfieldvalue('overlay_image')

    xlim = options.getfieldvalue('xlim', [min(md.mesh.x), max(md.mesh.x)])
    ylim = options.getfieldvalue('ylim', [min(md.mesh.y), max(md.mesh.y)])

    gtif = gdal.Open(image)
    trans = gtif.GetGeoTransform()
    xmin = trans[0]
    xmax = trans[0] + gtif.RasterXSize * trans[1]
    ymin = trans[3] + gtif.RasterYSize * trans[5]
    ymax = trans[3]
    # allow supplied image to have limits smaller than basemap or model limits
    x0 = max(min(xlim), xmin)
    x1 = min(max(xlim), xmax)
    y0 = max(min(ylim), ymin)
    y1 = min(max(ylim), ymax)
    inputname = 'temp.tif'
    os.system('gdal_translate-quiet - projwin ' + str(x0) + ' ' + str(y1) + ' ' + str(x1) + ' ' + str(y0) + ' ' + image + ' ' + inputname)

    gtif = gdal.Open(inputname)
    arr = gtif.ReadAsArray()
    #os.system('rm -rf . / temp.tif')

    if gtif.RasterCount >= 3:  # RGB array
        r = gtif.GetRasterBand(1).ReadAsArray()
        g = gtif.GetRasterBand(2).ReadAsArray()
        b = gtif.GetRasterBand(3).ReadAsArray()
        arr = 0.299 * r + 0.587 * g + 0.114 * b

    # normalize array
    arr = arr / np.float(np.max(arr.ravel()))
    arr = 1. - arr  # somehow the values got flipped

    if options.getfieldvalue('overlayhist', 0) == 1:
        ax = plt.gca()
        num = 2
        while True:
            if not plt.fignum_exists(num):
                break
            else:
                num += 1
        plt.figure(num)
        plt.hist(arr.flatten(), bins=256, range=(0., 1.))
        plt.title('histogram of overlay image, use for setting overlaylims')
        plt.show()
        plt.sca(ax) # return to original axes/figure

    # get parameters from cropped geotiff
    trans = gtif.GetGeoTransform()
    xmin = trans[0]
    xmax = trans[0] + gtif.RasterXSize * trans[1]
    ymin = trans[3] + gtif.RasterYSize * trans[5]
    ymax = trans[3]
    dx = trans[1]
    dy = trans[5]

    xarr = np.arange(xmin, xmax, dx)
    yarr = np.arange(ymin, ymax, -dy)  # - dy since origin = 'upper' (not sure how robust this is)
    xg, yg = np.meshgrid(xarr, yarr)
    overlaylims = options.getfieldvalue('overlaylims', [min(arr.ravel()), max(arr.ravel())])
    norm = mpl.colors.Normalize(vmin=overlaylims[0], vmax=overlaylims[1])

    pc = ax.pcolormesh(xg, yg, np.flipud(arr), cmap=mpl.cm.Greys, norm=norm)

    if options.exist('basemap'):
        # create coordinate grid in map projection units (for plotting)
        if md.mesh.epsg == 3413:
            hemisphere = 1
            st_lat = 70
            lon_0 = 45
        elif md.mesh.epsg == 3031:
            hemisphere = -1
            st_lat = 71
            lon_0 = 0
        else:
            hemisphere = eval(input('epsg code {} is not supported chose your hemisphere (1 for North, -1 for south)'.format(md.mesh.epsg)))

        lat, lon = xy2ll(xlim, ylim, hemisphere, lon_0, st_lat)
        extent = [np.diff(xlim)[0], np.diff(ylim)[0]]
        center = [lon[0] + np.diff(lon)[0] * 0.5, lat[0] + np.diff(lat)[0] * 0.5]
        m = Basemap(llcrnrlon=lon[0], llcrnrlat=lat[0], urcrnrlon=lon[1], urcrnrlat=lat[1],
                    lon_0=center[0], lat_0=center[1], width=extent[0], height=extent[1],
                    epsg=md.mesh.epsg, anchor='NW',
                    resolution='i', ax=ax)

        meridians = np.arange(np.floor(lon[0]), np.ceil(lon[1]), 1.)
        parallels = np.arange(np.floor(lat[0]), np.ceil(lat[1]), 1.)
        m.drawparallels(parallels, labels=[1, 0, 0, 0], ax=ax)  # labels = [left, right, top, bottom]
        m.drawmeridians(meridians, labels=[0, 0, 1, 0], ax=ax)
        m.drawcoastlines(ax=ax)
        m.drawmapboundary(ax=ax)

    #rasterization?
    if options.getfieldvalue('rasterized', 0):
        pc.set_rasterized(True)
