'use strict';

/**
 * Global Function: plotOverlay
 *
 * Plot a georeferenced image 
 *
 * Dependencies:
 *      - /ext/proj4-<version>
 *
 * See Also:
 *      - /js/plotmodel.js
 *      - /js/plotmanager.js
 */
function plot_overlay(md, data, options, canvas) {
    let elements        = null;
    let meshresults     = null; // results from call to processmesh
    let newX            = null; // scaled x-coordinates (array)
    let newY            = null; // scaled y-coordinates (array)
    let newZ            = null; // scaled z-coordinates (array)
    let position        = null; // x-, y-, and z-coordinates of vertex (object)
    let renderObject    = null; 
    let reproject       = false; // default; whether or not to reproject overlay
    let result          = null; // First, projection of radar overlay onto mesh; then, projection of overlay onto mesh
    let scale           = null; // coordinate scalars (array)
	let vertices        = null;
    let x               = null; // x-coordinates from results (array)
    let xlim            = null; // x-cooridinate limits (array)
    let y               = null; // y-coordinates from results (array)
    let ylim            = null; // y-cooridinate limits (array)
    let z               = null; // z-coordinates from results (array)
    
	// Process data and model
    meshresults = processmesh(md, [], options);
	x           = meshresults[0]; 
	y           = meshresults[1]; 
	z           = meshresults[2]; 
	elements    = meshresults[3];
	
	if (md.mesh.classname().startsWith('mesh2d')) {
		if (vesl.helpers.isNaN(md.geometry.surface) || (md.geometry.surface[0] !== undefined && vesl.helpers.isNaN(md.geometry.surface[0]))) {
			md.geometry.surface = NewArrayFill(x.length, 0);
			z = NewArrayFill(x.length, 0);
		} else {
			z = md.geometry.surface;
		}
	}
	
	// Compute coordinates and data range
	xlim = options.getfieldvalue('xlim', [ArrayMin(x), ArrayMax(x)]);
	ylim = options.getfieldvalue('ylim', [ArrayMin(y), ArrayMax(y)]);
	
	// Handle radaroverlay
	if (md.radaroverlay.outerx) {
		result      = Node.prototype.mergeVertices(x, y, z, elements, md.radaroverlay.outerx, md.radaroverlay.outery, md.radaroverlay.outerheight, md.radaroverlay.outerindex);
		x           = result.x;
		y           = result.y;
		z           = result.z;
		elements    = result.elements;
	}
	
	// Handle heightscale
	renderObject = options.getfieldvalue('render', {});
	
	if ('overlay' in renderObject && renderObject.overlay.reproject !== undefined) {
		reproject = renderObject.overlay.reproject;
	}
	
	if (reproject) {
		//NOTE: This is being hardcoded for now, as we want to display Josh's greenland-ice-retreat model on a globe.
		// If mesh3dprisims do not need to be displayed on a globe, add an option to the view/render objects to specify the reprojection.
		// Currently, we are also hardcoding the reprojection from EPSG:3413 (WGS 84 / NSIDC Sea Ice Polar Stereographic North) to EPSG:4326 (WGS 84, standard 2d spherical projection))
		// projection definitions are taken from http://epsg.io/XXXXX.js, where XXXXX is the projection number.
		// Can be dynamically taken and retrieved via md.mesh.epsg, but again, hardcoded for efficiency.
		proj4.defs("EPSG:3413","+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"); // http://epsg.io/3413.js
		//Don't need to redefine EPSG:4326, as it is defined by default - however, this is still included just in case.
		// proj4.defs("EPSG:4326","+proj=longlat +datum=WGS84 +no_defs"); // http://epsg.io/4326.js
		//Create empty copies of XYZ to avoid rescaling md.geomtery.surface in the case that z = md.geomtery.surface
		newX = NewArrayFill(x.length, 0);
		newY = NewArrayFill(x.length, 0);
		newZ = NewArrayFill(x.length, 0);
		
		for (let i = 0; i < x.length; ++i) {
			result      = proj4('EPSG:3413','EPSG:4326', [x[i], y[i]]);
			position    = vesl.geo.geographicToCartesian(vesl.EARTH_RADIUS + z[i], result[1], result[0], true);
			newX[i]     = position.x;
			newY[i]     = position.y;
			newZ[i]     = position.z;
		}
		
		x           = newX;
		y           = newY;
		z           = newZ;
		vertices    = Node.prototype.scaleVertices(md, x, y, z, elements, options.getfieldvalue('heightscale', 1), options.getfieldvalue('maskregion',{'enabled':false, 'reproject':reproject}), true);
		scale       = [1, 1, 1];
    } else {
		vertices    = Node.prototype.scaleVertices(md, x, y, z, elements, options.getfieldvalue('heightscale', 1), options.getfieldvalue('maskregion',{'enabled':false, 'reproject':reproject}), true);
		scale       = [1, 1, 1];
	}
	
	// Compute coordinates and data range
	xlim = options.getfieldvalue('xlim', [ArrayMin(x), ArrayMax(x)]);
	ylim = options.getfieldvalue('ylim', [ArrayMin(y), ArrayMax(y)]);
	
	// Compute GL variables
	let texture         = new THREE.TextureLoader().load(options.getfieldvalue('overlay_image'));
	let groundEnabled   = false;
	let name            = 'overlay';
	
	if ('ground' in renderObject) {
    	groundEnabled = renderObject.ground.enabled;
    }
    
	if ('overlay' in renderObject) {
    	name = renderObject.overlay.name;
    }
    
	let shaderName = groundEnabled ? 'ground' : 'overlay';
	
	let node = new Node(
		'canvas', canvas,
		'options', options,
		'md', md,
		'name', name,
		'shaderName', shaderName,
		'opacity', options.getfieldvalue('outeropacity', 1.0),
		'cullFace', THREE.DoubleSide,
		'texture', texture,
		'scale', scale
	);
	canvas.overlayNode = node;

	let xRange      = xlim[1] - xlim[0];
	let yRange      = ylim[1] - ylim[0];
	let coordArray  = [new Array(x.length), new Array(x.length)];
	
	// Generate mesh
	if (reproject || md.mesh.classname().startsWith('mesh3d')) {
		let magnitude   = 0;
		let xyz         = null;
		
		for (let i = 0; i < x.length; ++i) {
			xyz                 = vec3.fromValues(vertices[0][i], vertices[1][i], vertices[2][i]);
			magnitude           = vec3.length(xyz);
		
			coordArray[0][i]    = Math.atan2(xyz[1], xyz[0]) / (2 * Math.PI) + 0.5;
			coordArray[1][i]    = Math.asin(xyz[2] / magnitude) / Math.PI + 0.5;
		}
	} else {
		for (let i = 0; i < x.length; ++i) {
			coordArray[0][i] = (vertices[0][i] - xlim[0]) / xRange;
			coordArray[1][i] = (vertices[1][i] - ylim[0]) / yRange;
		}
	}
	
	node.patch('Faces', elements, 'Vertices', vertices, 'FaceVertexCData', coordArray, 'FaceColor', 'interp');
}
