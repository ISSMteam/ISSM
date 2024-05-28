'use strict';

function applyoptions(md, data, options, canvas) {
	//APPLYOPTIONS - apply colobar, text, cloud, and expdisp options to current plot
	//
	//   Usage:
	//      applyoptions(md, data, options)
	//
	//   See also: PLOTMODEL, PARSE_OPTIONS

	//{{{ colorbar
	let gl = canvas.gl;
	
	if (options.getfieldvalue('colorbar', false)) {
		//{{{ Create colorbar labels
		let cAxis = options.getfieldvalue('caxis');
		let labels = [];
		let divisions = options.getfieldvalue('colorbarnticks', 6);
		let cAxisDelta = cAxis[1] - cAxis[0];
		let precision = options.getfieldvalue('colorbarprecision', 3);
		let format = options.getfieldvalue('colorbarformat', 'f').toLowerCase();
		if (options.getfieldvalue('log','off') !== 'off') {
			for (let i = divisions; i >= 0; i--) {
				let scale = (Math.log10(cAxis[1]) - Math.log10(cAxis[0])) / Math.log10(options.getfieldvalue('log', 10));
				if (format === 'f') {
					labels[i] = (Math.pow(options.getfieldvalue('log', 10), Math.log10(cAxis[0]) / Math.log10(options.getfieldvalue('log', 10)) + scale * (divisions - i) / divisions)).toFixed(precision);
				} else if (format === 'e') {
					labels[i] = (Math.pow(options.getfieldvalue('log', 10), Math.log10(cAxis[0]) / Math.log10(options.getfieldvalue('log', 10)) + scale * (divisions - i) / divisions)).toPrecision(precision);
				} else {
					labels[i] = (Math.pow(options.getfieldvalue('log', 10), Math.log10(cAxis[0]) / Math.log10(options.getfieldvalue('log', 10)) + scale * (divisions - i) / divisions)).toFixed(precision);
				}
			}
		} else {
			for (let i = divisions; i >= 0; i--) {
				if (format === 'f') {
					labels[i] = (cAxisDelta * (divisions - i) / divisions + cAxis[0]).toFixed(precision);
				} else if (format === 'e') {
					labels[i] = (cAxisDelta * (divisions - i) / divisions + cAxis[0]).toPrecision(precision);
				} else {
					labels[i] = (cAxisDelta * (divisions - i) / divisions + cAxis[0]).toFixed(precision);
				}
			}
		} //}}}
		//{{{ Initialize colorbar canvas
		let cCanvasId = options.getfieldvalue('colorbarid', options.getfieldvalue('canvasid') + ('-colorbar-canvas'));
		let cCanvasIdBase = cCanvasId.substring(0, cCanvasId.lastIndexOf('-canvas'));
		let cCanvas = document.getElementById(cCanvasId);
		let cWidth = cCanvas.width * options.getfieldvalue('colorbarwidth', 1);
		let cHeight = cCanvas.height * options.getfieldvalue('colorbarheight', 1);
		let cContext = cCanvas.getContext('2d');
		let cMap = options.getfieldvalue('colormap', 'amp');
		
		// If value of cMap is of type array, assume we have a custom colorbar; otherwise, look up colormap by name from global letiable "colorbars"
		let colorbar = null;
		
		if (vesl.arrays.isArray(cMap)) {
			colorbar = cMap;
		} else {
    		for (let colormap in colorbars) {
        		if (colormap === cMap) {
			        colorbar = colorbars[cMap];
			        break;
                }
            }
            
            if (colorbar === null) {
                for (let colormap in cmoceanColormaps) {
            		if (colormap === cMap) {
    			        colorbar = cmoceanColormaps[cMap];
    			        break;
                    }
                }
            }
		}
		
		let gradient = cContext.createLinearGradient(0, 0, 0, cHeight);
		//}}}
		//{{{ Draw colorbar gradient
		// TODO: Allow for passing the opacity in as a fourth value of each array of a colormap?
		let applyOpacityToColorbar  = options.getfieldvalue('applyOpacityToColorbar', false);
		let color                   = null;
		let background              = options.getfieldvalue('colorbarBackground', null);
		let offset                  = 1 / (colorbar.length - 1) / 2;
		let opacity                 = options.getfieldvalue('opacity', 1.0);
		let position                = null;
		let scaling                 = 1 - 2 * offset;
		
		if (!applyOpacityToColorbar) {
    		opacity = 1.0;
		}
		
		if (background !== null) {
    		background = [background[0] * 255, background[1] * 255, background[2] * 255];
    		$('#' + cCanvasId).css('background', 'rgba(' + background.toString() + ', 1.0');
		}
		
		for (let i = 0; i < colorbar.length; i++) {
			color = colorbar[colorbar.length - i - 1];
			color = [Math.round(color[0] * 255), Math.round(color[1] * 255), Math.round(color[2] * 255)];
			position = (i / (colorbar.length - 1) * scaling) + offset;
			gradient.addColorStop(position, 'rgba(' + color.toString() + ', ' + opacity + ')');
		}
		
		cContext.clearRect(0, 0, cWidth, cHeight);
		cContext.beginPath();
		cContext.fillStyle = gradient;
		cContext.fillRect(0, 0, cWidth, cHeight);
		//}}}
		//{{{ Draw colorbar border
		cContext.beginPath();
		cContext.lineWidth = '1';
		cContext.strokeStyle=options.getfieldvalue('colorbarfontcolor','black');
		cContext.rect(0, 0, cWidth, cHeight);
		cContext.stroke();
		//}}}
		//{{{ Draw colorbar labels
		let cLabelsId = cCanvasIdBase + '-labels';
		let cLabels = $('#' + cLabelsId);
		let cLabelString = '';
		let x, y;
		cLabels.empty();
		for (let i = 0; i <= divisions; i++) {
			y = (i + 0.5) / (divisions + 1) * cHeight;
			x = 0.2 * cWidth;
			cLabelString += '<li><span>' + labels[i] + '</span></li>';
			cContext.beginPath();
			cContext.moveTo(0, y);
			cContext.lineTo(x, y);
			cContext.moveTo(cWidth - x, y);
			cContext.lineTo(cWidth, y);
			cContext.stroke();
		}
		cLabels.append(cLabelString);
		//}}}
		//{{{ Draw colorbar title
		let cTitleId = cCanvasIdBase + '-heading';
		let cTitle = $('#' + cTitleId);
		if (options.exist('colorbartitle')) { cTitle.html(options.getfieldvalue('colorbartitle')); }
		//}}}
		//{{{ Setup texture/alpha canvases
		let $canvas 	= $(canvas);
		let tCanvasId 	= options.getfieldvalue('texturecanvasid', 'texturecanvas');
		let aCanvasId 	= options.getfieldvalue('alphacanvasid', 'alphacanvas');
		let tCanvas 	= document.getElementById(tCanvasId);
		let aCanvas 	= document.getElementById(aCanvasId);
		
		if (tCanvas == null) {
			$('<canvas id="' + tCanvasId + '" width="256" height="256" style="display: none;"></canvas>').insertAfter($canvas);
			tCanvas = document.getElementById(tCanvasId);
		}
		
		if (aCanvas == null) {
			$('<canvas id="' + aCanvasId + '" width="256" height="256" style="display: none;"></canvas>').insertAfter($canvas);
			aCanvas = document.getElementById(aCanvasId);
		}
	
		//Set up canvas drawing contexes and gradients.
		let tContext = tCanvas.getContext('2d');
		let aContext = aCanvas.getContext('2d');
		let tGradient = tContext.createLinearGradient(0, 0, 0, 256);
		let aGradient = aContext.createLinearGradient(0, 0, 0, 256);
		
		//Determine where in gradient to start unit mesh transparency
		let maskAlphaEnabled = options.getfieldvalue('maskAlphaEnabled', false);
		let maskAlphaTolerance = options.getfieldvalue('maskAlphaTolerance', 0.1);
		let maskAlphaValue = options.getfieldvalue('maskAlphaValue', 1.1);
		let maskAlphaUseColor = options.getfieldvalue('maskAlphaUseColor', false);
		let maskAlphaColor = options.getfieldvalue('maskAlphaColor', 'rgba(0.0, 0.0, 255, 1.0)');
		let alphaValue = (maskAlphaValue - cAxis[0]) / cAxisDelta;
		
		//Apply transparency to alpha map that enables alpha to be read from texture, and to actual texture alpha.
		for (let i = 0; i < colorbar.length; i++) {
			color = colorbar[colorbar.length - i - 1];
			color = [Math.round(color[0] * 255), Math.round(color[1] * 255), Math.round(color[2] * 255)];
			let colorStop = i / (colorbar.length - 1);
			if (maskAlphaEnabled && (colorStop > 1 - alphaValue || colorStop == colorbar.length - 1)) {
				if (maskAlphaUseColor) {
					tGradient.addColorStop(colorStop, maskAlphaColor);
					aGradient.addColorStop(colorStop, 'rgb(255, 255, 255)');
				} else {
					tGradient.addColorStop(colorStop, 'rgba(' + color.toString() + ', 0.0)');
					aGradient.addColorStop(colorStop, 'rgb(0, 0, 0)');
				}
			} else {
				tGradient.addColorStop(colorStop, 'rgba(' + color.toString() + ', 1.0)');
				aGradient.addColorStop(colorStop, 'rgb(255, 255, 255)');
			}
		}
		
		//Draw gradients to canvaes.
		tContext.fillStyle = tGradient;
		aContext.fillStyle = aGradient;
		tContext.fillRect(0, 0, 256, 256);
		aContext.fillRect(0, 0, 256, 256);
		
		//Allow for special texture colors, drawing each color in equal width vertical rectangles. The last rectanglar section is reserved for the colormap.
		if (options.exist('maskregion')) {
			let maskObject = options.getfieldvalue('maskregion',{'enabled':false});
			if (maskObject.enabled && !vesl.helpers.isEmptyOrUndefined(maskObject.colors)) {
				let x = 0;
				let sections = Object.keys(maskObject.colors).length + 1;
				let size = 256;
				let width = Math.floor(1 / sections * size);
				for (let color in maskObject.colors) {
					tContext.fillStyle = maskObject.colors[color];
					tContext.fillRect(x++ * width, 0, width, size);
				}
			}
		}
		
		//Read canvases as images, and load as textures in Three.js
		let tURL            = tCanvas.toDataURL();
		let aURL            = aCanvas.toDataURL();
		let textureMap      = new THREE.TextureLoader().load(tURL);
		let alphaMap        = new THREE.TextureLoader().load(aURL);
		let unitOptions 	= options.getfieldvalue('unitOptions', {'name' : 'unit'});
		let unitName 		= unitOptions.name;
		let unitSceneNode 	= canvas.unitNodes[unitName].sceneNode;
		unitSceneNode.material.map          = textureMap;
		unitSceneNode.material.emissiveMap  = textureMap;
		unitSceneNode.material.color        = new THREE.Color(0xffffff);
		unitSceneNode.material.needsUpdate  = true;
		
		//Only apply alpha map if enabled.
		if (maskAlphaEnabled) {
			unitSceneNode.material.alphaMap = alphaMap;
		}
	} //}}}
	//}}}
	//{{{ Data marker
	// TODO: Default parameters are already being handled by /js/vesl/DataMarker.js::constructor, so perhaps we should be initializing this elsewhere and/or by some other manner	
	if (vesl.helpers.isEmptyOrUndefined(canvas.dataMarker)) { // Only define data marker once
		let dataMarker 			= {};
		let dataMarkerOptions 	= options.getfieldvalue('dataMarker', {});
		
		dataMarker = {
			enabled: defaultFor(dataMarkerOptions.enabled, false)
		};
		
		// Only initialize data marker object if we have enabled data markers for this canvas
		if (dataMarker.enabled) {
			dataMarker.object = new vesl.DataMarker(
				{
					canvas 			: canvas,
					hasMarker 		: defaultFor(dataMarkerOptions.hasMarker, true),
					markerImgPath 	: defaultFor(dataMarkerOptions.markerImgPath, '/canvas/data-markers/data-marker.svg'),
					hasTooltip		: defaultFor(dataMarkerOptions.hasTooltip, false),
					tooltipFormat	: defaultFor(dataMarkerOptions.tooltipFormat, ''),
					tooltipFields	: defaultFor(dataMarkerOptions.tooltipFields, null),
					width			: defaultFor(dataMarkerOptions.width, 32),
					height			: defaultFor(dataMarkerOptions.height, 32)
				}
			);
		}
		
		canvas.dataMarker = dataMarker;
	}
	//}}}
	//contours
	if (options.exist('contourlevels')) {
		plot_contour(md,data,options,canvas);
	}
} //}}}
function drawGroundingLines(md, canvas, options, renderObject, lines, colors) { //{{{
	let renderObjects = options.getfieldvalue('render',{});
	let state = canvas.state;
	let scene = state.scene;
	
	let group = scene.getObjectByName('groundingLines');
	if (group !== undefined) {
		scene.remove(group); //Remove old group if already exists
	}
	group = new THREE.Group();
	group.name = 'groundingLines';
	
	//Plot multiple grounding lines, each consisting of multiple polygons. 
	for (let i = 0; i < lines.length; i++) {
		let groundingLine = lines[i]; //of type Proxy (not object or array), thus must iterate without length attribute
		let color = processColor(colors[i]);
		let x = [];
		let y = [];
		let z = [];
		//In order to show polygons correctly, must convert from line strip to line segments
		for (let j in groundingLine) {
			let polygon = groundingLine[j]
			let lineStripX = polygon['x'];
			let lineStripY = polygon['y'];
			let lineStripZ = polygon['z'];
			let lineSegmentsX = [lineStripX[0]];
			let lineSegmentsY = [lineStripY[0]];
			let lineSegmentsZ = [lineStripZ[0]];
			//Must push same coordinates as end of previous segment and beginning of next segment
			for (let k = 1; k < lineStripX.length - 1; k++) {
				lineSegmentsX.push(lineStripX[k], lineStripX[k]);
				lineSegmentsY.push(lineStripY[k], lineStripY[k]);
				lineSegmentsZ.push(lineStripZ[k], lineStripZ[k]);
			}
			//Cap off last coordinate
			lineSegmentsX.push(lineStripX[lineStripX.length - 1]);
			lineSegmentsY.push(lineStripY[lineStripY.length - 1]);
			lineSegmentsZ.push(lineStripZ[lineStripZ.length - 1]);
			//Add polygon coordinates to existing groundingLine
			x = x.concat(lineSegmentsX);
			y = y.concat(lineSegmentsY);
			if (renderObject.followsBed) {
				z = z.concat(lineSegmentsZ);
			} else {
				z = z.concat(NewArrayFill(lineSegmentsX.length, 10 * i));
			}
		}
		z = ArrayScale(z, options.getfieldvalue('heightscale', 1));
		let vertices = [x, y, z];
		//console.log('ArrayMin(x):', ArrayMin(x), 'ArrayMax(x):', ArrayMax(x), 'ArrayMin(y):', ArrayMin(y), 'ArrayMax(y):', ArrayMax(y));
			
		let node = new Node(
			'canvas', canvas,
			'options', options,
			'renderObject', renderObject,
			'name', 'groundingLine_' + String(i),
			'shaderName', 'lines',
			//'dashed', dashed,
			'depthTest', false,
			'diffuseColor', color,
			'lineWidth', renderObject.lineWidth,
			//'lineWidth', options.getfieldvalue('linewidth', 3),
			'scale', [renderObject.scale, renderObject.scale, renderObject.scale]
		);
		node.patch('Vertices', vertices, 'group', group, 'ViewReset', false);
	}
	scene.add(group);
	return group;
} //}}}
function drawText(state, options, textLabels) { //{{{
	let camera = state.camera;
	let overlayCanvas = state.overlayCanvas;
	let ctx = overlayCanvas.getContext('2d');
	let widthHalf = overlayCanvas.width / 2
	let heightHalf = overlayCanvas.height / 2;
	
	for (let i = 0; i < textLabels.length; i++) {
		let textLabel = textLabels[i];
		let position = textLabel.position.clone();
		let text = textLabel.text;
		
		//Project world coordinates to screenspace coordinates
		position.project(camera);
		let x = (position.x * widthHalf) + widthHalf;
		let y = -(position.y * heightHalf) + heightHalf;
		
		ctx.font = 'bold ' + String(options.getfieldvalue('colorbarfontsize', 22))+'px Arial Black, sans-serif';
		ctx.fillStyle = options.getfieldvalue('colorbarfontcolor','black');
		ctx.strokeStyle = 'white';
		ctx.textAlign = 'center';
		ctx.textBaseline = 'middle';
		ctx.fillText(text, x, y);
		ctx.strokeText(text, x, y);
	}
} //}}}
function processColor(color) { //{{{
	//Update the diffuse color with an RGB color name or vec4 containing r, g, b, and opacity values from 0.0 to 1.0
	if (typeof color === 'string') {
		color = new THREE.Color(color);
	} else if (Array.isArray(color)) {
		color = new THREE.Color().fromArray(color);
	}
	return color;
} //}}}
