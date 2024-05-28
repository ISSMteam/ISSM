# Test Name: Sea-Level-Partitions
#
# TODO:
# - Save boundaries by name to some data structure with tag so that they can be
#   copied to basins that use identical shapefiles and projections (will need
#   to check if the cost of the additional structure, checks, and copying are
#   greater than the cost of projecting).
#
from averaging import averaging
from BamgTriangulate import BamgTriangulate
from basin import *
from epsg2proj import epsg2proj
from find_point import find_point
from GetAreas3DTria import GetAreas3DTria
from InterpFromMesh2d import InterpFromMesh2d
from InterpFromMeshToMesh2d import InterpFromMeshToMesh2d
from intersect import intersect
from lovenumbers import *
from materials import materials
from model import *
from plotmodel import plotmodel
from sealevelmodel import *
from solve import solve


testagainst2002 = 0

# Data paths {{{
shppath = '../Data/shp/'
gshhsshapefile = shppath + 'GSHHS_c_L1-NightlyRun.shp'
#}}}

#create sealevel model to hold our information
sl = sealevelmodel()

#Create basins using boundaries from shapefile
#some projections we'll rely on #{{{
proj4326 = epsg2proj(4326)
proj3031 = epsg2proj(3031)
#}}}
#HemisphereWest #{{{
sl.addbasin(
    basin('continent', 'hemispherewest', 'name', 'hemispherewest', 'proj', laea(0, -90), 'boundaries', [  # Peru projection 3587
        boundary('shppath', shppath, 'shpfilename', 'HemisphereSplit', 'proj', proj4326, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'NorthAntarctica', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneBrunt', 'proj', proj3031, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'RonneEastSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneFront', 'proj', proj3031, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'RonneWestSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'WestAntarctica2', 'proj', proj3031, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'SouthAntarctica', 'proj', proj3031)]
    )
)
#}}}
#Ross: {{{
sl.addbasin(
    basin('continent', 'antarctica', 'name', 'ross', 'proj', proj3031, 'boundaries', [
        boundary('shppath', shppath, 'shpfilename', 'SouthAntarctica', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossIceShelf', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossWestFront', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossFront', 'proj', proj3031, 'orientation', 'reverse')]
    )
)
#}}}
#HemisphereEast: {{{
sl.addbasin(
    basin('continent', 'hemisphereeast', 'name', 'hemisphereeast', 'proj', laea(0, +90), 'boundaries', [  #Australian projection lat, long
        boundary('shppath', shppath, 'shpfilename', 'HemisphereSplit', 'proj', proj4326),
        boundary('shppath', shppath, 'shpfilename', 'SouthAntarctica', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossFront', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossWestFront', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'EastAntarctica2', 'proj', proj3031, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'NorthAntarctica', 'proj', proj3031)]
    )
)
#}}}
#Antarctica excluding Ronne: {{{
sl.addbasin(
    basin('continent', 'antarctica', 'name', 'antarctica-grounded', 'proj', proj3031, 'boundaries', [
        boundary('shppath', shppath, 'shpfilename', 'NorthAntarctica', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'EastAntarctica2', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossWestFront', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RossIceShelf', 'proj', proj3031, 'orientation', 'reverse'),
        boundary('shppath', shppath, 'shpfilename', 'SouthAntarctica', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'WestAntarctica2', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneWestSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneIceShelf', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneEastSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneBrunt', 'proj', proj3031)]
    )
)
#}}}
#Ronne: {{{
sl.addbasin(
    basin('continent', 'antarctica', 'name', 'ronne', 'proj', proj3031, 'boundaries', [
        boundary('shppath', shppath, 'shpfilename', 'RonneWestSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneIceShelf', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneEastSummit', 'proj', proj3031),
        boundary('shppath', shppath, 'shpfilename', 'RonneFront', 'proj', proj3031, 'orientation', 'reverse')]
    )
)
#}}}

# Meshing
# Go through basins and mesh #{{{
# Meshing parameters: {{{
hmin = 500
hmax = 700
hmin = hmin * 1000
hmax = hmax * 1000
tolerance = 100  # tolerance of 100m on Earth position when merging 3d meshes
threshold = 5
defaultoptions = ['KeepVertices', 0,
                  'MaxCornerAngle', 0.0000000001,
                  'NoBoundaryRefinement', 1]
alreadyloaded = 0
#}}}
for ind in sl.basinindx('basin', 'all'):
    bas = sl.basins[ind]
    print('Meshing basin {}'.format(bas.name))

    # Recover basin domain
    domain = bas.contour()

    # Recover coastline inside basin, using GSHHS_c_L1. It's a lat/long file, hence epsg 4326
    coastline = bas.shapefilecrop('shapefile', gshhsshapefile, 'epsgshapefile', 4326, 'threshold', threshold)

    # Mesh
    md = bamg(model(), 'domain', domain, 'subdomains', coastline, 'hmin', hmin, 'hmax', hmax, *defaultoptions)  # NOTE: Unpacking defaultoptions with '*'

    # Miscellaneous
    md.mesh.proj = bas.proj
    md.miscellaneous.name = bas.name

    # Recover mask where we have land
    md.private.bamg.landmask = (md.private.bamg['mesh'].Triangles[:, 3] >= 1).astype(int)

    # Vertex connectivity
    md.mesh.vertexconnectivity = NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)

    # Add model to sl icecaps
    sl.addicecap(md)
#}}}

# Parameterization
# Parameterize ice sheets #{{{
for ind in sl.basinindx('continent', ['antarctica']):
    print('Parameterizing basin {}'.format(sl.icecaps[ind].miscellaneous.name))

    md = sl.icecaps[ind]
    bas = sl.basins[ind]

    # Masks #{{{
    # Ice levelset from domain outlines
    md.mask.ice_levelset = -np.ones((md.mesh.numberofvertices, ))

    if bas.isnameany('antarctica-grounded'):
        md.mask.ocean_levelset = np.ones((md.mesh.numberofvertices, ))

    if bas.isnameany('ronne', 'ross'):
        md.mask.ocean_levelset = -np.ones((md.mesh.numberofvertices, ))
    #}}}

    # Lat/long #{{{
    md.mesh.long, md.mesh.lat = gdaltransform(md.mesh.x, md.mesh.y, md.mesh.proj, 'EPSG:4326')
    #}}}

    # Geometry #{{{
    if bas.iscontinentany('antarctica'):
        di = md.materials.rho_ice / md.materials.rho_water

        print('      reading bedrock')
        md.geometry.bed = -np.ones((md.mesh.numberofvertices, ))
        md.geometry.base = md.geometry.bed
        md.geometry.thickness = 1000 * np.ones((md.mesh.numberofvertices, ))
        md.geometry.surface = md.geometry.bed + md.geometry.thickness
    #}}}

    # SLC #{{{
    if bas.iscontinentany('antarctica'):
        if testagainst2002:
            # TODO: Check if the following works as expected: 'pos' is empty, so nothing is assigned to 'md.solidearth.surfaceload.icethicknesschange[pos]'
            md.masstransport.spcthickness = np.zeros((md.mesh.numberofvertices, ))
            # Antarctica
            late = np.sum(md.mesh.lat[md.mesh.elements - 1], axis=1) / 3
            longe = np.sum(md.mesh.long[md.mesh.elements - 1], axis=1) / 3
            pos = np.where(late < -85)[0]
            ratio = 0.225314032985172 / 0.193045366574523
            md.masstransport.spcthickness[md.mesh.elements[pos]] = md.masstransport.spcthickness[md.mesh.elements[pos]] - 100 * ratio
        else:
            delH = np.loadtxt('../Data/AIS_delH_trend.txt')
            longAIS = delH[:, 0]
            latAIS = delH[:, 1]
            delHAIS = delH[:, 2]
            index = BamgTriangulate(longAIS, latAIS)
            late = md.mesh.lat
            longe = md.mesh.long + 360
            pos = np.where(longe > 360)[0]
            longe[pos] = longe[pos] - 360
            delHAIS = InterpFromMesh2d(index, longAIS, latAIS, delHAIS, longe, late) # NOTE: Compare to corresponding output under MATLAB to understand why we offset triangle indices by 1 (only caught because Triangle.cpp was producing triangles with negative areas)
            northpole = find_point(md.mesh.long, md.mesh.lat, 0, 90)
            delHAIS[northpole] = 0
            md.masstransport.spcthickness = delHAIS / 100

        md.initialization.sealevel = np.zeros((md.mesh.numberofvertices, ))

        md.dsl.global_average_thermosteric_sea_level = np.zeros((2, 1))
        md.dsl.sea_surface_height_above_geoid = np.zeros((md.mesh.numberofvertices + 1, 1))
        md.dsl.sea_water_pressure_at_sea_floor = np.zeros((md.mesh.numberofvertices + 1, 1))
    #}}}

    # Material properties #{{{
    md.materials = materials('hydro')
    #}}}

    # Diverse #{{{
    md.miscellaneous.name = bas.name
    #}}}

    sl.icecaps[ind] = md
#}}}

# Parameterize continents #{{{
for ind in sl.basinindx('continent', ['hemisphereeast', 'hemispherewest']):
    print('Masks for basin {}'.format(sl.icecaps[ind].miscellaneous.name))
    md = sl.icecaps[ind]
    bas = sl.basins[ind]

    # Recover lat, long
    md.mesh.long, md.mesh.lat = gdaltransform(md.mesh.x, md.mesh.y, md.mesh.proj, 'EPSG:4326')

    # Mask #{{{
    # Figure out mask from initial mesh: deal with land and ocean masks (land includes grounded ice). #{{{
    # First, transform land element mask into vertex-driven one
    land = md.private.bamg.landmask
    land_mask = -np.ones((md.mesh.numberofvertices, ))

    landels = np.nonzero(land)[0]
    land_mask[md.mesh.elements[landels, :] - 1] = 1 # NOTE: Need to offset each element of each row of md.mesh.elements[landels, :] by one

    # Go through edges of each land element
    connectedels = md.mesh.elementconnectivity[landels, :] - 1 # NOTE: Need to offset each element of each row of md.mesh.elementconnectivity[landels, :] by one
    connectedisonocean = np.logical_not(land[connectedels]).astype(int)
    sumconnectedisonocean = np.sum(connectedisonocean, axis=1)

    # Figure out which land elements are connected to the ocean
    landelsconocean = landels[np.nonzero(sumconnectedisonocean)[0]]

    ind1 = np.hstack((md.mesh.elements[landelsconocean, 0],
                      md.mesh.elements[landelsconocean, 1],
                      md.mesh.elements[landelsconocean, 2]))
    ind2 = np.hstack((md.mesh.elements[landelsconocean, 1],
                      md.mesh.elements[landelsconocean, 2],
                      md.mesh.elements[landelsconocean, 0]))

    # Edge ind1 and ind2
    for i in range(len(ind1)):
        els1 = md.mesh.vertexconnectivity[ind1[i] - 1, 0:md.mesh.vertexconnectivity[ind1[i] - 1, -1]]
        els2 = md.mesh.vertexconnectivity[ind2[i] - 1, 0:md.mesh.vertexconnectivity[ind2[i] - 1, -1]]
        els = intersect(els1, els2)[0] # NOTE: Throwing away second- and third- position values returned from call

        if len(np.nonzero(land[els - 1])[0]) == 1:
            # This edge is on the beach, 0 to the edge
            land_mask[ind1[i] - 1] = 0
            land_mask[ind2[i] - 1] = 0

    md.mask.ocean_levelset = land_mask
    md.mask.ice_levelset = np.ones((md.mesh.numberofvertices, )) # If there are glaciers, we'll adjust

    if testagainst2002:
        #{{{
        # Greenland
        pos = np.where(np.logical_and.reduce((md.mesh.lat > 70, md.mesh.lat < 80, md.mesh.long > -60, md.mesh.long < -30)))[0]
        md.mask.ice_levelset[pos] = -1
        #}}}
    #}}}

    # SLC loading/calibration #{{{
    md.masstransport.spcthickness = np.zeros((md.mesh.numberofvertices, ))

    if testagainst2002:
        #{{{
        # Greenland
        late = np.sum(md.mesh.lat[md.mesh.elements - 1], axis=1) / 3
        longe = np.sum(md.mesh.long[md.mesh.elements - 1], axis=1) / 3
        pos = np.where(np.logical_and.reduce((late > 70, late < 80, longe > -60, longe < -30)))[0]
        ratio = .3823 / .262344
        md.masstransport.spcthickness[md.mesh.elements[pos]] = md.masstransport.spcthickness[md.mesh.elements[pos]] - 100 * ratio

        # Correct mask
        md.mask.ice_levelset[md.mesh.elements[pos, :] - 1] = -1
        #}}}
    else:
        delH = np.loadtxt('../Data/GIS_delH_trend.txt')
        longGIS = delH[:, 0]
        latGIS = delH[:, 1]
        delHGIS = delH[:, 2]
        index = BamgTriangulate(longGIS, latGIS)
        late = md.mesh.lat
        longe = md.mesh.long + 360
        pos = np.where(longe > 360)[0]
        longe[pos] = longe[pos] - 360
        delHGIS = InterpFromMeshToMesh2d(index, longGIS, latGIS, delHGIS, longe, late)

        delH = np.loadtxt('../Data/GLA_delH_trend_15regions.txt')
        longGLA = delH[:, 0]
        latGLA = delH[:, 1]
        delHGLA = np.sum(delH[:, 2:], axis=1)
        index = BamgTriangulate(longGLA, latGLA)
        late = md.mesh.lat
        longe = md.mesh.long + 360
        pos = np.where(longe > 360)[0]
        longe[pos] = longe[pos] - 360
        delHGLA = InterpFromMeshToMesh2d(index, longGLA, latGLA, delHGLA, longe, late)

        # NOTE: For some reason, cannot apply pos to multiple arrays in a
        #       singlelike we might do in MATLAB. Instead, we iterate over
        #       elements of pos.
        #
        pos = np.nonzero(delHGIS)[0]
        for p in pos:
            md.masstransport.spcthickness[p] = md.masstransport.spcthickness[p] - delHGIS[p] / 100
        pos = np.nonzero(delHGLA)[0]
        for p in pos:
            md.masstransport.spcthickness[p] = md.masstransport.spcthickness[p] - delHGLA[p] / 100

        # Adjust mask accordingly
        pos = np.nonzero(md.masstransport.spcthickness)[0]
        md.mask.ice_levelset[pos] = -1
        md.mask.ocean_levelset[pos] = 1

    md.initialization.sealevel = np.zeros((md.mesh.numberofvertices, ))

    md.dsl.global_average_thermosteric_sea_level = np.zeros((2, 1))
    md.dsl.sea_surface_height_above_geoid = np.zeros((md.mesh.numberofvertices + 1, 1))
    md.dsl.sea_water_pressure_at_sea_floor = np.zeros((md.mesh.numberofvertices + 1, 1))
    #}}}

    # Geometry #{{{
    di = md.materials.rho_ice / md.materials.rho_water
    md.geometry.bed = -np.ones((md.mesh.numberofvertices, ))
    md.geometry.base = md.geometry.bed
    md.geometry.thickness = 1000 * np.ones((md.mesh.numberofvertices, ))
    md.geometry.surface = md.geometry.bed + md.geometry.thickness
    #}}}
    # Materials #{{{
    md.materials = materials('hydro')
    #}}}
    sl.icecaps[ind] = md

# Assemble Earth in 3D #{{{

# Parameters
plotting = 0
tolerance = 100
loneedgesdetect = 0

# Create Earth model by concatenating all the icecaps in 3D
sl.caticecaps('tolerance', tolerance, 'loneedgesdetect', loneedgesdetect)

# Figure out how each icecap's mesh connects to the larger Earth mesh
sl.intersections('force', 1)

# Figure out connectivity
print('Mesh connectivity')
sl.earth.mesh.vertexconnectivity = NodeConnectivity(sl.earth.mesh.elements, sl.earth.mesh.numberofvertices)

# Areas
print('Mesh nodal areas')
sl.earth.mesh.area = averaging(sl.earth, GetAreas3DTria(sl.earth.mesh.elements, sl.earth.mesh.x, sl.earth.mesh.y, sl.earth.mesh.z), 4)

# Transfer a list of fields from each icecap and continent back to Earth
sl.transfer('mask.ice_levelset')
sl.transfer('mask.ocean_levelset')
sl.transfer('geometry.bed')
sl.transfer('geometry.surface')
sl.transfer('geometry.thickness')
sl.transfer('geometry.base')
sl.transfer('mesh.lat')
sl.transfer('mesh.long')
sl.transfer('masstransport.spcthickness') #
sl.transfer('initialization.sealevel')
sl.transfer('dsl.sea_surface_height_above_geoid')
sl.transfer('dsl.sea_water_pressure_at_sea_floor')

# Radius
sl.earth.mesh.r = (sl.earth.mesh.x ** 2 + sl.earth.mesh.y ** 2 + sl.earth.mesh.z ** 2) ** 0.5 # NOTE: math.sqrt cannot be applied element-wise to a list/numpy.array

# Check on the mesh transitions #{{{
plotting = 0
if plotting:
    flags = np.ones((sl.earth.mesh.numberofelements, ))
    for i in range(len(sl.eltransitions)):
        flags[sl.eltransitions[i]] = i
    # TODO: Reconfigure the following in the process of bringing plotting online
    plotmodel(sl.earth, 'data', flags, 'shading', 'faceted', 'coastline', 'on', 'coast_color', 'g')
#}}}

#}}}

# Solve Sea-level equation on Earth only #{{{
md = sl.earth # we don't do computations on ice sheets or land

#Materials
md.materials = materials('hydro')

# Elastic loading from love numbers
md.solidearth.lovenumbers = lovenumbers('maxdeg', 100)
md.solidearth.settings.ocean_area_scaling = 0

# Miscellaneous
md.miscellaneous.name = 'test2004'

# New stuff
md.dsl.global_average_thermosteric_sea_level = np.array([[(1.1 + .38)], [0]]) # steric + water storage AR5

# Solutuion parameters
md.solidearth.settings.reltol = np.nan
md.solidearth.settings.abstol = 1e-3
md.solidearth.settings.sealevelloading = 1
md.solidearth.settings.isgrd = 1
md.solidearth.settings.ocean_area_scaling = 1
md.solidearth.settings.grdmodel = 1
md.timestepping.time_step = 1

# Physics
md.transient.issmb = 0
md.transient.isstressbalance = 0
md.transient.isthermal = 0
md.transient.ismasstransport = 1
md.transient.isslc = 1

# Initializations
md.basalforcings.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices,))
md.basalforcings.floatingice_melting_rate = np.zeros((md.mesh.numberofvertices,))
md.initialization.vx = np.zeros((md.mesh.numberofvertices,))
md.initialization.vy = np.zeros((md.mesh.numberofvertices,))
md.initialization.sealevel = np.zeros((md.mesh.numberofvertices,))
md.initialization.bottompressure = np.zeros((md.mesh.numberofvertices,))
md.initialization.dsl = np.zeros((md.mesh.numberofvertices,))
md.initialization.str = 0
md.smb.mass_balance = np.zeros((md.mesh.numberofvertices,))

# Max number of iterations reverted back to 10 (i.e. the original default value)
md.solidearth.settings.maxiter = 10

# Eustatic run:
md.solidearth.settings.selfattraction = 0
md.solidearth.settings.elastic = 0
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md.solidearth.requested_outputs = ['default',
                                   'DeltaIceThickness',
                                   'Sealevel',
                                   'Bed',
                                   'SealevelBarystaticIceMask',
                                   'SealevelBarystaticOceanMask']
md = solve(md, 'Transient')
Seustatic = md.results.TransientSolution.Sealevel

# Eustatic + selfattraction run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 0
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md = solve(md, 'Transient')
Sselfattraction = md.results.TransientSolution.Sealevel

# Eustatic + selfattraction + elastic run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 0
md.solidearth.settings.viscous = 0
md = solve(md, 'Transient')
Selastic = md.results.TransientSolution.Sealevel

# Eustatic + selfattraction + elastic + rotation run
md.solidearth.settings.selfattraction = 1
md.solidearth.settings.elastic = 1
md.solidearth.settings.rotation = 1
md.solidearth.settings.viscous = 0
md = solve(md, 'Transient')
Srotation = md.results.TransientSolution.Sealevel

#Fields and tolerances to track changes
field_names = ['Eustatic', 'Rigid', 'Elastic', 'Rotation']
field_tolerances = [1e-13, 1e-13, 1e-13, 1e-13]
field_values = [Seustatic, Sselfattraction, Selastic, Srotation]
