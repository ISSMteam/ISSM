The Matlab scripts in the kml directory may be divided into three areas:
    1. kml class definitions
    2. kml utilities
    3. kml drivers

Each of the three will be discussed below.

I.  KML Class Definitions

All the kml is defined using the reference:
http://code.google.com/apis/kml/documentation/kmlreference.html

Many of the kml objects described in the first figure of the kml reference are defined as Matlab objects, with most or all of the kml properties.  The classes consist of one base class and several levels of derived classes (super and sub, respectively, in Matlab terminology).  The hierarchy of these classes are as follows.

kml_object.m
    kml_feature.m
        kml_placemark.m
        kml_container.m
            kml_folder.m
            kml_document.m
    kml_geometry.m
        kml_point.m
        kml_linestring.m
        kml_linearring.m
        kml_polygon.m
        kml_multigeometry.m
    kml_styleselector.m
        kml_style.m
    kml_substyle.m
        kml_colorstyle.m
            kml_linestyle.m
            kml_polystyle.m

For each class, the methods consist of a constructor, a display method, a fieldnames method (since objects can not use the Matlab "fieldnames" function), a set method, and a write method.  All of the methods, as well as the properties, inherit the base classes where appropriate rather than repeating any functionality.  The Matlab "help" command can be used to list the documentation for any of the classes.

II.  KML Utilities

The kml utilities take an ISSM model (and optionally results) and write it into a kml object, typically a kml folder, in various ways.  The kml object could then be written into a file, by itself or with others, and imported into Google Earth.  Note that kml polygons can only be one color, not interpolated.

The following kml utilities have been written.  The kml_mesh_elem and kml_partitions are probably the most useful, because the rest are intermediate results.  The Matlab "help" command can again be used to list the documentation.

kml_mesh_elem.m      - write a kml folder with each ISSM element as a kml polygon (color-coded by results, if provided)
kml_part_flagedges.m - write a kml folder with each segment between two ISSM partitions as a kml linestring
kml_unsh_edges.m     - write a kml folder with each unshared segment of an ISSM partition as a kml linestring
kml_part_elems.m     - write a kml folder with all the elements of each ISSM partition as kml polygons, noting that elements are repeated for each partition in which they have nodes (color-coded by results, if provided)
kml_part_edges.m     - write a kml folder with all the edges of each ISSM partition as a kml linestring (color-coded by results, if provided)
kml_partitions.m     - write a kml folder with each ISSM partition as a kml polygon (color-coded by results, if provided)

In order the write any and all kml objects that have been constructed to a file for import into Google Earth, the following function has been written.

kml_file_write.m     - write a kml file of the specified kml objects (may open and/or close file)

III.  KML Drivers

There is one kml driver, which may be used as a model for custom drivers.  It takes an ISSM model, converts the data from node to element (if necessary), and writes some kml headers and style templates.  In addition, it constructs six kml folders for the first six kml utilities above, then writes those to the file and closes the file.

kml_mesh_write.m     - write a kml file of the ISSM model (color-coded by results, if provided)

IV.  Other Utilities

There are some other utilities that are used in the construction of topological tables for the kml writing.

kmlnodeconnectivity.m   - create a node connectivity table (nnodes x mxepg+1)
edgeadjacency.m      - create an edge adjacency array (elems x edges)
edgeperimeter.m      - create an edge perimeter (edgeper x 2) and element perimeter (edgeper x 1) list

