function project2d(md3d,value,layer) {
    // PROJECT2D - returns the value of a field for a given layer of the mesh
    // 
    //    returns the value of a vector for a given layer from extruded mesh onto the 2d mesh 
    //    used to do the extrusion. This function is used to compare values between different
    //    layers of a 3d mesh.
    // 
    //    Usage:
    //       projection_value=project2d(md3d,value,layer)
    // 
    //    Example:
    //       vel2=project2d(md3d,md3d.initialization.vel,2);
    //       returns the velocity of the second layer (1 is the base)

    // some checks on list of arguments
    if (arguments.length !== 3) {
        console.error('project2d error message');
    }

    if (md3d.mesh.domaintype() !== '3D') {
        console.error("wrong model type ... should be ''3d''");
    }

    if (layer<1 || layer>md3d.mesh.numberoflayers) {
        console.error(['layer must be between 1 and ' + num2str(md3d.mesh.numberoflayers)]);
    }

    // Return the projection value
    var temp = [];
    if (value.length === md3d.mesh.numberofvertices) {
        for (var i = (layer-1)*md3d.mesh.numberofvertices2d; i <= layer*md3d.mesh.numberofvertices2d; ++i) {
            temp.push(value[i]);
        }
    } else if (value.length === md3d.mesh.numberofvertices+1) {
        for (var i = (layer-1)*md3d.mesh.numberofvertices2d; i <= layer*md3d.mesh.numberofvertices2d; ++i) {
            temp.push(value[i]);
        }
        temp.push(value[value.length-1]);
    } else {
        for (var i = (layer-1)*md3d.mesh.numberofelements; i <= layer*md3d.mesh.numberofelements2d; ++i) {
            temp.push(value[i]);
        }
    }

    return temp;
}
