function DepthAverage(md,vector) {
    // DEPTHAVERAGE - computes depth average of 3d vector using the trapezoidal rule, and returns the value on 2d mesh. 
    // 
    //    Usage:
    //       vector_average=DepthAverage(md,vector);
    // 
    //    Example:
    //       vel_bar=DepthAverage(md,md.initialization.vel);

    // check that the model given in input is 3d
    if (md.mesh.elementtype() !== 'Penta') {
        console.error('DepthAverage error message: the model given in input must be 3d');
    }

    // nods data
    if (vector.length === md.mesh.numberofvertices) {
        var vector_average=zeros(md.mesh.numberofvertices2d,1);

        for (var i = 1; i < md.mesh.numberoflayers-1; ++i) {
            vector_average = vector_average.map(function(x) {
                return x + (project2d(md, vector, i) + project2d(md,vector,i+1))/2;
            }).map(function(y) {
                return y * (project2d(md, md.mesh.z, i+1) - project2d(md, md.mesh.z, i));
            });
        }

        vector_average = vector_average.map(function(z) {
            return z / project2d(md, md.geometry.thickness, 1);
        });

        return vector_average;
    }
    // element data
    else if (vector.length === md.mesh.numberofelements) {
        var vector_average=zeros(md.mesh.numberofelements2d,1);
        for (var i = 1; i < md.mesh.numberoflayers-1; ++i) {
            vector_average = vector_average.map(function(x) {
                return x + project2d(md, vector, i);
            }).map(function(y) {
                return y * (project2d(md, md.mesh.z, i+1) - project2d(md, md.mesh.z, i));
            });
        }

        vector_average = vector_average.map(function(z) {
            return z / project2d(md, md.geometry.thickness, 1);
        });

        return vector_average;
    } else {
        console.error('vector size not supported yet');
    }
}
