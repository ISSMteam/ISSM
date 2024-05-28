from plotoptions import plotoptions
from checkplotoptions import checkplotoptions
from model import model
import numpy as np
import math
from writejsfile import writejsfile


def export_gl(md, *varargin):
    class ResultObj(object):
        def __getattr__(self, attr):
            return self.__dict__.get(attr)

    print('getting options')
    templist = plotoptions(varargin)
    optionslist = templist.list
    options = optionslist[1]
    options = checkplotoptions(md, options)
    #print (templist, options)
    #templist contains options 0 - 3. Use in the future to rework.

    #Setup unique directory in present dir:
    print('setting directory')
    directory = optionslist[0].getfieldvalue('directory')
    databasename = optionslist[0].getfieldvalue('database')

    #scaling factor:
    print('setting scaling factor')
    scaling_factor = optionslist[0].getfieldvalue('scaling_factor')

    #Deal with title:
    print('setting title')
    if optionslist[0].exist('title'):
        title = optionslist[0].getfieldvalue('title')
    else:
        title = ''

    #initialize model:
    print('initializing model')
    model.title = title
    model.initialZoomFactor = options.getfieldvalue('zoom', -.25)

    #Deal with contour {{{
    print('getting contour')
    print(md.mesh.segments)
    segmenets0 = [s - 1 for s in md.mesh.segments[:, 0]]
    segmenets1 = [s - 1 for s in md.mesh.segments[:, 1]]

    contour_lat1 = md.mesh.lat.take(segmenets0)
    contour_lat2 = md.mesh.lat.take(segmenets1)
    contour_long1 = md.mesh.long.take(segmenets0)
    contour_long2 = md.mesh.long.take(segmenets1)
    contour_surface1 = md.geometry.surface.take(segmenets0)
    contour_surface2 = md.geometry.surface.take(segmenets1)

    R1 = 6371000 * np.ones(len(contour_surface1)) + scaling_factor * contour_surface1
    R2 = 6371000 * np.ones(len(contour_surface2)) + scaling_factor * contour_surface2

    model.contourx1 = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.cos(math.radians(int)), R1, contour_lat1, contour_long1))
    model.contoury1 = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.sin(math.radians(int)), R1, contour_lat1, contour_long1))
    model.contourz1 = list(map(lambda r, lat: r * math.sin(math.radians(lat)), R1, contour_lat1))

    model.contourx2 = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.cos(math.radians(int)), R2, contour_lat2, contour_long2))
    model.contoury2 = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.sin(math.radians(int)), R2, contour_lat2, contour_long2))
    model.contourz2 = list(map(lambda r, lat: r * math.sin(math.radians(lat)), R2, contour_lat2))

    # }}}
    #Deal with mesh and results {{{
    print('getting mesh')
    surface = md.geometry.surface.flatten()
    numberofelements = md.mesh.numberofelements
    numberofvertices = md.mesh.numberofvertices
    R = 6371000 * np.ones(len(md.mesh.lat)) + scaling_factor * surface

    x = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.cos(math.radians(int)), R, md.mesh.lat, md.mesh.long))
    y = list(map(lambda r, lat, int: r * math.cos(math.radians(lat)) * math.sin(math.radians(int)), R, md.mesh.lat, md.mesh.long))
    z = list(map(lambda r, lat: r * math.sin(math.radians(lat)), R, md.mesh.lat))

    #Deal with triangulation:
    print('getting triangulation')
    model.index = md.mesh.elements
    model.x = x
    model.y = y
    model.z = z
    model.surface = surface

    results = []
    print(optionslist)
    #Deal with data:
    print('getting data')
    for i in range(0, len(optionslist)):
        options = optionslist[i]
        options = checkplotoptions(md, options)
        data = options.getfieldvalue('data').flatten()
        results.append(ResultObj())
        results[i].data = data
        results[i].caxis = options.getfieldvalue('caxis', [min(data), max(data)])

        label = options.getfieldvalue('label', '')
        if label == '':
            #create generic label:
            label = ['data', str(i)]
        results[i].label = label

        shortlabel = options.getfieldvalue('shortlabel', '')
        if shortlabel == '':
            #create generic short label:
            shortlabel = ['data', str(i)]
        results[i].shortlabel = shortlabel

        if type(data[2]) != np.float64:
            time_range = options.getfieldvalue('time_range', [0, 100])
            results[i].time_range = time_range

        unit = options.getfieldvalue('unit', '')
        if unit == '':
            #create generic unit:
            unit = 'SI'
        results[i].unit = unit
    model.results = results

    #Write model to javascript database file:
    print('writing to file')
    writejsfile(directory + databasename + '.js', model, databasename)
    # }}}
