/*
	DIRECTIVES
*/

/* globals jQuery */

// Execute script in strict mode
'use strict';

/******************************************************************************/

/**
 * expread - Read a string containing the contents of an .exp file and return an object
 *
 *	This function takes as input a string containing the contents of an .exp file read in with FileReader.readAsText(),
 *	and an instance of class Callout, which is used for error reporting to the user.
 *
 *	The function builds and outputs an array of objects, each containing,
 *		the file name,
 *		the number of nodes,
 *		the density, 
 *		and a boolean representing whether or not the domain is closed
 *	for each of the profiles represented in the input file.
 *
 *	Usage:
 *		let contours = expread(file, callout)
 *
 *	Example (assumes HTML5, jQuery, and class vesl.Callout):
 *      HTML:
 *			<div id="loadDomainCalloutId" class="callout">
 *				<h4 class="callout-header"></h4>
 *				<div class="callout-content"></div>
 * 			</div>
 *	
 *			<input id="fileInput" type="file">
 *		
 *		Javascript:
 *			let loadDomainCallout = new vesl.Callout('loadDomainCalloutId');
 *
 *			$('#fileInput').change(function(event) {
 *  			let contours 	= {};
 *	 			let file 		= event.target.files[0];
 *				let fileReader 	= new FileReader();
 *							
 *				fileReader.onload = function(event) {
 *					contours = expread(event.target.result, loadDomainCallout);
 *                  // Now if contours is not empty or undefined, do something with it
 *				};
 *						
 *				fileReader.readAsText(file);
 *			});
 *
 *	See also /js/Callout.js
 */
function expread(file, callout) {//{{{
	return new Promise(function(resolve, reject) {
		/*
			Constants
		*/
		//{{{
		let CALLOUT_ERROR_HEADER = 'Oh no!';
		//}}}
		
		
		/*
			Variables
		*/
		//{{{
		let contour 	= {};
		let contours 	= [];
		let count 		= 0;
		let lineTokens 	= [];
		let lines 		= [];
		let linesIndex	= 0;
		//}}}
		
		
		// Split file contents on either Linux or Windows line breaks
		lines = file.split(/[\r\n]+/g);
		
		// Loop over the number of profiles
		while (lines[linesIndex + 1] !== undefined) { // May need to convert this comparison if parsing of EOF as undefined changes in the future
			// Update number of profiles
			contour 			= {};
			contours[count++] 	= contour; 
			
			/*
				Get file name
			*/
			//{{{
			lineTokens = lines[linesIndex++].split(/[\s\t]+/);
			
			if (!(lineTokens.length === 2 && lineTokens[0] === '##' && lineTokens[1].slice(0, 5) === 'Name:')) {
				callout.set(CALLOUT_ERROR_HEADER, 'File name line of profile ' + count + ' is not in a valid format.');
				callout.show();
				reject('File name line of profile ' + count + ' is not in a valid format.');
			}
			
			if (lineTokens[1].length > 5) {
				contour['name'] = lineTokens[1].slice(5);
			} else {
				contour['name'] = '';
			}
			//}}}
			
			
			/*
				Get icon
			*/
			//{{{
			lineTokens = lines[linesIndex++].split(/[\s\t]+/);
			
			if (!(lineTokens.length === 2 && lineTokens[0] === '##' && lineTokens[1].slice(0, 5) === 'Icon:')) {
				callout.set(CALLOUT_ERROR_HEADER, 'Icon line of profile ' + count + ' is not in a valid format.');
				callout.show();
				reject('Icon line of profile ' + count + ' is not in a valid format.');
			}
			//}}}
			
			
			/*
				Get info
			*/
			//{{{
			lineTokens = lines[linesIndex++].split(/[\s\t]+/);
			
			if (!(lineTokens.length === 4 && lineTokens[0] === '#' && lineTokens[1] === 'Points')) {
				callout.set(CALLOUT_ERROR_HEADER, 'First info line of profile ' + count + ' is not in a valid format.');
				callout.show();
				reject('First info line of profile ' + count + ' is not in a valid format.');
			}
			//}}}
	
			
			// Get number of nodes and density
			//{{{
			lineTokens = lines[linesIndex++].split(/[\s\t]+/);
			
			contour['nods'] 	= parseInt(lineTokens[0]);
			contour['density'] 	= parseInt(lineTokens[1]);
			//}}}
			
			
			// Get info
			//{{{
			lineTokens = lines[linesIndex++].split(/[\s\t]+/);
			
			if (!(lineTokens.join('') === '#XposYpos')) {
				callout.set(CALLOUT_ERROR_HEADER, 'Second info line of profile ' + count + ' is not in a valid format.');
				callout.show();
				reject('Second info line of profile ' + count + ' is not in a valid format.');
			}
			//}}}
			
			
			// Get coordinates
			//{{{
			contour['x'] = [];
			contour['y'] = [];
			
			for (let i = 0; i < contour['nods']; ++i) {
				lineTokens = lines[linesIndex++].split(/[\s\t]+/); // Note that after last iteration of this loop, lines[linesIndex] points to blank line following the most recently parsed profile		
				contour['x'][i] = parseFloat(lineTokens[0]);
				contour['y'][i]	= parseFloat(lineTokens[1]);
			}
			//}}}
			
			
			// Check if closed
			if (contour['nods'] > 1 &&
				contour['x'][-1] === contour['x'][0] &&
				contour['y'][-1] === contour['y'][0]) {
				contour['closed'] = true;	
			} else {
				contour['closed'] = false;
			}
		}
		
		resolve(contour);
	});
}//}}}
