function plotdoc() { //{{{
	//PLOTDOC - plot documentation
	//
	//   Usage:
	//      plotdoc()
	
	//TODO: rename innermask/outermask/maskzero and combine

	console.log('   Plot usage: plotmodel(model,varargin)');
	console.log('   Options: ');
	console.log('       "canvasid": canvas id');
	console.log('       "data" : what we want to plot');
	console.log('       	Available values for "data" are: ');
	console.log('       		- any field of the model structure. ex: plot(md,"data","vel"), or plot(md,"data",md.initialization.vel)');
	console.log('       		- "mesh": draw mesh using trisurf');
	console.log('       		- "quiver": quiver plot');
	console.log('       "backgroundcolor": plot background color. (default "lightcyan", ex: "green","blue")');
	console.log('       "brush": specify brush options (default {"strength":0.075,"falloff":0.5})');
	console.log('       	"enabled": toggle brush (default false, ex: true)');
	console.log('       	"strength": value that brush will change data points by (ex: 0.075)');
	console.log('       	"falloff": multiplier that brush will decrease strength by for each successive point away from brush center (ex: 0.5)');
	console.log('       "clouds": specify brush options (default {"strength":0.075,"falloff":0.5})');
	console.log('       	"enabled": toggle clouds (default false, ex: true)');
	console.log('       	"height": height to spawn clouds at (ex: 7500)');
	console.log('       	"quantity": quantity of clouds to spawn (ex: 10)');
	console.log('       "caxis": modify colorbar range (array of type [a, b] where b>=a)');
	console.log('       "colorbar": add colorbar (default "off", ex: "on", "off")');
	console.log('       "colorbarid": colorbar canvas id (string)');
	console.log('       "colorbartitle": colorbar title (string)');	
	console.log('       "colorbarnticks": number of colorbar ticks (default 6, ex: 2, 10)');	
	console.log('       "colorbarprecision": colorbar label digit precision (default 3, ex: 1, 4)');	
	console.log('       "colorbarinnerlabels": choose whether labels are inside colorbar (default "off", ex: "on", "off")');
	console.log('       "colorbarfontsize": specify colorbar font size (default 1, ex: 14, 22)');
	console.log('       "colorbarfontcolor": specify colorbar font color (default "black", ex: "green","blue")');
	console.log('       "colorbarwidth": multiplier (default 1) to the default width colorbar');
	console.log('       "colorbarheight": multiplier (default 1) to the default height colorbar');
	console.log('       "colormap": same as standard matlab option (default "jet", ex: "hsv","cool","spring","gray","Ala","Rignot",...)');
	console.log('       "controlsensitivity": sensitivty of view/zoom changes as a percentage of default (default 1, ex: 0.5, 2.75)');
	console.log('       "dataMarker": object cotaining data marker parameters. See webgl.js for defaults. (ex: {"enabled":true,"format":["<div id="sim-plot"></div>"],"labels":["thickness","velocity","value"],"animated":true})');
	console.log('       	"enabled": toggle data marker displays (default true, ex: false)');
	console.log('       	"hasMarker": whether or not the data marker has a marker (default true, ex: false)');
	console.log('       	"markerImgPath": path to image to use for marker (default "/canvas/data-markers/data-marker.svg")');
	console.log('       	"hasTooltip": whether or not the data marker has a tooltip (default false, ex: true)');
	console.log('       	"tooltipFormat": C-style format string to be used when on calls to vesl.DataMarker.tooltipContent to format passed data points (default "", ex: ["X: %.2e<br>Y: %.2e<br>Z: %.2e])');	
	console.log('       	"width": width, in pixels, of the data marker (default 32)');	
	console.log('       	"height": height, in pixels, of the data marker (default 32)');
	console.log('       "displayview": print view value to console (default "off", ex: "on", "off")');
	console.log('       "displayzoom": print zoom value to console (default "off", ex: "on", "off")');
	console.log('       "edgecolor": same as standard matlab option EdgeColor (default "black", ex: color name: "blue" or RGB array: [0.5, 0.2, 0.8])');
	console.log('       "heightscale": scaling factor to accentuate height. (default 1, ex: 0.5, 100)');
	console.log('       "linewidth": line width for mesh, quiver, and contour plots, currently limited by WebGL to 1. (default 1, ex: 2, 5)');
	console.log('       "log": value of log (default 10, ex: 2, Math.E)');
	console.log('       "mask": list of flags of size numberofnodes or numberofelements. Only "true" values are plotted ');
	console.log('       "movies": object cotaining transient plot animation options (ex: {"fps":4,"loop":true})');
	console.log('       "maskzeros": object cotaining transient plot animation options (ex: "enabled":true,"color":[1.0, 1.0, 1.0, 1.0],"tolerance":1e-3,"zeroValue":0.5})');
	console.log('       	"enabled": toggle maskzeros (default false, ex: true)');
	console.log('       	"color": RGBA color value array with ranges 0.0 to 1.0 (ex: [1.0, 1.0, 1.0, 1.0])');
	console.log('       	"tolerance": values within this tolerance of the zeroValue will be masked. (default: 1e-3, ex: 2.5e-2)');
	console.log('       	"zeroValue": the percentage value with range 0.0, to 1.0 of the caxis value around which the data will be masked with the color. (default: 0.5, ex: 0, 1.0, 0.75)');
	console.log('       "innermask*": Special mask that colors all parts of a data mesh below a height a certain color. provide innermaskheight and innermaskcolor options also (default "off", ex: "on", "off")');
	console.log('       "outermask*": Special mask that colors all parts of a overlay mesh below a height a certain color. provide outermaskheight and outermaskcolor options also (default "off", ex: "on", "off")');
	console.log('       "overlay": overlay a radar amplitude image behind (default "off", ex: "on", "off")');
	console.log('       "overlay_image": path to overlay image (default "", ex: "./img/radar.png")');
	console.log('       "quiver": add quiver plot overlay for velocities. (default "off", ex: "on", "off")');
	console.log('       "scaling": scaling factor used by quiver plots. Default is 0.4');
	console.log('       "opacity": transparency coefficient 0.0 to 1.0, the lower, the more transparent. (default 1.0, ex: 0.5, 0.25)');
	console.log('       "render": a object containing a list of default object to render. (default {}, ex: {"sky", "space"})');
	console.log('       	"sky": render the atmosphere. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle sky (default false, ex: true)');
	console.log('       	"space": render space. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle space (default false, ex: true)');
	console.log('       	"coastlines": render coastlines. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle coastlines (default false, ex: true)');
	console.log('       		"scale": scale coastlines factor (default 1.0, ex: 1.004)');
	console.log('       		"x": x coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"y": y coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"z": z coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       	"city": render city. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle city (default false, ex: true)');
	console.log('       		"size": radius of city sphere, in meters (default 1.0, ex: 150000)');
	console.log('       		"color": color of city sphere (ex: "magenta")');
	console.log('       		"x": x coordinate of city. (ex: 0.0)');
	console.log('       		"y": y coordinate of city. (ex: 0.0)');
	console.log('       		"z": z coordinate of city. (ex: 6356700.0)');
	console.log('       	"cities": render cities. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle cities (default false, ex: true)');
	console.log('       		"size": radius of cities spheres, in meters (default 1.0, ex: 80000)');
	console.log('       		"color": color of cities spheres (ex: "darkviolet")');
	console.log('       		"x": x coordinate array of cities. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"y": y coordinate array of cities. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"z": z coordinate array of cities. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       	"graticule": render graticule. (ex: {"enabled":true})');
	console.log('       		"enabled": toggle graticule (default false, ex: true)');
	console.log('       		"scale": scale graticule factor (default 1.0, ex: 1.004)');
	console.log('       		"x": x coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"y": y coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       		"z": z coordinate array. (ex: [0.0, 10.0, -25.0,...])');
	console.log('       "view": object cotaining view parameters. See webgl.js for defaults. (ex: {"position":[0.0,0.0,0.0],"rotation":[0.0,0.0,0.0],"zoom":1.0,"zoomLimits":[0.01,100.0],"azimuthLimits":[-180,180.0],"elevationLimits":[-90,90.0],"panningEnabled":false,"twod":false})');
	console.log('       	"position": camera position (ex: [0.0,0.0,0.0])');
	console.log('       	"rotation": camera rotation (ex: [0.0,0.0,0.0])');
	console.log('       	"zoom": initial camera zoom as a percentage of default (default 1, ex: 1.5, 0.01)');
	console.log('       	"zoomLimits": zoom view limits (ex: [0.05, 10])');
	console.log('      	 	"azimuthLimits": zoom view limits (ex: [0.05, 10])');
	console.log('       	"elevationLimits": zoom view limits (ex: [0.05, 10])');
	console.log('       	"panningEnabled": controls panning with shift + drag mouse or pan gestures (default: false, ex: true)');
	console.log('       	"twod": controls twod orthographic view (default: false, ex: true)');
	console.log('       "xlim": x coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "ylim": y coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "zlim": z coordinates to fit inside camera (ex: [0, 500])');
	console.log('       "transient_field_data": array of data objects (ex: [[0.0,1.0, 2.5, 12.0...],[0.0,1.0, 2.5, 12.0...],...])');
	console.log('       "textlabels": plot text labels rendered in 3d space, using an array of text/coordinate pairs (ex: [{"pos":[0.0,0.0,0.0],"text":"origin"}])');
	
	console.log('  ');
	console.log('   Examples:');
	console.log('       plotmodel(md,"data","vel","data","mesh","view#2",3,"colorbar#all","on")');
	console.log('       plotmodel(md,"data",md.geomtery.surface)');
} //}}}
