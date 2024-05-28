// Execute script in strict mode
'use strict';
/**
 * generic - Class that allows for sending requests to server and handling response or errors
 *
 *	Usage:
 *		generic = new generic();
 *
 *	Todo:
 *		Convert to ES6 class
 */
function generic() {
	// Retrieve options and apply to properties 
	// {{{
	let args 	= Array.prototype.slice.call(arguments);
	let options = new pairoptions(args.slice(0, args.length));
	
	this.url			= options.getfieldvalue('url', '');
	this.np				= options.getfieldvalue('np', 3);
	this.codeversion	= options.getfieldvalue('codeversion', 20486);
	this.codepath		= options.getfieldvalue('codepath', 'issmdir/bin');
	this.executionpath	= options.getfieldvalue('executionpath', 'issmdir/execution');
	//}}}
	
	
	// Methods
	//{{{
	this.disp = function() { //{{{
		console.log(sprintf('   generic class echo:'));
		console.log(sprintf('    url: %s', this.url));
		console.log(sprintf('    np: %i', this.np));
		console.log(sprintf('    codepath: %s', this.codepath));
		console.log(sprintf('    executionpath: %s', this.executionpath));
	}; //}}}
	
	this.classname = function() { //{{{
		return 'generic';
	}; //}}}
	
	this.checkconsistency = function(md, solution, analyses) { //{{{
		if (cluster.np < 1) {
			md.checkmessage('Number of processors should be at least 1!');
		}
		if (isNaN(cluster.np)) {
			md.checkmessage('Number of processors NaN!');
		}
	}; //}}}
	
	this.BuildQueueScript = function(cluster, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota) { //{{{
		//write queuing script 
		//what is the executable being called? 
		executable 	= 'issm.exe';
		
		fid 		= fopen(modelname + '.queue','w');
		fprintf(fid, '#!%s\n', cluster.shell);
		fprintf(fid, 'mpiexec -np %i %s/%s %s %s %s 2> %s.errlog >%s.outlog ', cluster.np, cluster.codepath, executable, EnumToString(solution), cluster.executionpath + '/' + dirname, modelname, modelname, modelname);					
		fclose(fid);
	}; //}}}
	
	this.UploadAndRun = function(md, fid, toolkitsstring, solutionstring, name, runtimename, successCallback, errorCallback, solveButtonId, callout, withProgressBar) { //{{{
		/* Local constants */
		let PROGRESS_BAR_ID 				= 'solve-progress-bar';
		let PROGRESS_BAR_TEXT_PERCENTAGE_ID = 'progress-bar-text-percentage';		
		
		/* Local variables */
		let hasCallout 					= false;
		let isProgressBarLoaded			= false;
		let progressBar					= {};
		let progressBarTextPercentage 	= {};
		let request 					= {};
		let solveButton 				= $(solveButtonId);
		let solveButtonText 			= !vesl.helpers.isEmptyOrUndefined(solveButton) ? solveButton.text() : ''; // Save initial solve button text
		
		/* Local functions */
		// NOTE: After conversion of generic to ES6 class, these should be class methods
		function loadProgressBar() {
			callout.setContent('\
				<div class="progress-bar-wrapper">\
					<progress id="' + PROGRESS_BAR_ID + '" value="0" max="100"></progress>\
					<div class="progress-bar-text">\
						<span id="' + PROGRESS_BAR_TEXT_PERCENTAGE_ID + '">0</span>\
						<span>%</span>\
					</div>\
				</div>\
			');
			
			progressBar 				= $('#' + PROGRESS_BAR_ID);
			progressBarTextPercentage 	= $('#' + PROGRESS_BAR_TEXT_PERCENTAGE_ID);
			isProgressBarLoaded 		= true;
		}
		
		function setProgressBar(progress) {
			progressBar.val(progress);
			progressBarTextPercentage.text(progress);
		}
		
		// Check certain arguments
		hasCallout = !vesl.helpers.isEmptyOrUndefined(callout);		
		
		// Check that we have a connection
		if (!navigator.onLine) {
			console.log('Error: no connection!');
			if (hasCallout) {
				callout.set('No connection!', '');
			} else {
				solveButton.text('No connection!').prop('disabled', false);
			}
			errorCallback();
			return;
		}
		
		// Create request
		request = new XMLHttpRequest();
		request.open('POST', this.url, true);
		
		// Set request properties
		request.position 	 = 0; // Keep track of current parsing position in repsonseText
		request.timeout 	= 1000 * 60 * 10; // in milliseconds; NOTE: Under current implementation, this should match 'Timeout' option in Apache conf (/etc/apache2/apache2.conf)
		request.resultBegin = false;
		
		request.ontimeout = function(event) { //{{{
			console.log('Error: timeout!');
			
			if (hasCallout) {
				callout.set(vesl.ERROR_HEADER_GENERAL, '<p>Request timeout! ' + vesl.ERROR_CONTENT_TRY_AGAIN + '</p>');
			} else {
				solveButton.text('Timeout!').prop('disabled', false);
			}
			
			errorCallback();
		}; //}}}
		
		request.onerror = function(event) { //{{{
			console.log('Error: could not run!');
			
			if (hasCallout) {
				callout.set(vesl.ERROR_HEADER_GENERAL, '<p>Something went wrong. ' + vesl.ERROR_CONTENT_TRY_AGAIN + '</p>');
			} else {
				solveButton.text('Could not run!').prop('disabled', false);
			}
			
			errorCallback();
		}; //}}}
		
		request.upload.onprogress = function(event) { //{{{
			let progress = (event.loaded / event.total * 100).toFixed(0);
			
			if (hasCallout) {
				callout.setHeader('Sending request...');
				if (withProgressBar) {
					if (!isProgressBarLoaded) {
						loadProgressBar();
					}
					setProgressBar(progress);
				} else {
					callout.setContent('<p>' + progress + '%</p>');
				}
			} else {
				solveButton.text('Sending: ' + progress + '%');
			}
        }; //}}}
        
		request.onprogress = function(event) { //{{{
			/* Local variables */
			let progress 	= 0;

			// Receive updates by parsing message length as a 32-bit hex string of form 0x*09ABCDEF))
			let startIndex 	= request.position;
			let endIndex 	= request.position + 10;
			
			if (request.responseText.length >= endIndex) { // Ensure entire hex string is loaded
				let chunkSize = parseInt(request.responseText.slice(startIndex, endIndex));
				
				startIndex 	= endIndex;
				endIndex 	= startIndex + chunkSize;
				
				if (chunkSize >= 1024) { // Arbitrary maximium size of message (Must be below minimium size of model results)
					progress = ((request.responseText.length - request.position) / chunkSize * 100).toFixed(0);
					if (hasCallout) {
						callout.setHeader('Downloading result...');
						if (withProgressBar) {
							setProgressBar(progress);
						} else {
							callout.setContent('<p>' + progress + '%</p>');
						}
					} else {
						solveButton.text('Downloading: ' + progress + '%').prop('disabled', true);
					}
				} else if (request.responseText.length >= endIndex) { // Ensure entire chunk is loaded
					if (request.resultBegin) {
						return;
					}
					chunk = request.responseText.slice(startIndex, endIndex)
					if (chunk == 'RESULT_BEGIN') {
						progress = 100;
						request.resultBegin = true;
					} else {
						progress = parseInt(chunk);
					}
					if (hasCallout) {
						callout.setHeader('Computing...');
						if (withProgressBar) {
							setProgressBar(progress);
						} else {
							callout.setContent('<p>' + progress + '%</p>');
						}
					} else {
						solveButton.text('Computing: ' + progress + '%').prop('disabled', true);
					}
					request.position = endIndex;
				}
			}
		}; //}}}
		
		request.onload = function(event) { //{{{
			/* Local variables */
			//{{{
			let buffer 				= {};
// 			let bufferInflated		= {};
			let responseText 		= '';
			//}}}
			
			// TODO: Get context to this.str2ab or otherwise declare it in order to avoid duplication
			function str2ab(str) { //{{{
				let buffer = new Uint8Array(str.length);
				for (let i = 0, strLen = str.length; i < strLen; ++i) {
					buffer[i] = str.charCodeAt(i);
				}
				return buffer;
			} //}}}
			
			try {
				//Scan through buffer until progress updates stop and 'RESULT_BEGIN' string is received
				let startIndex 	= request.position;
				let endIndex 	= request.position + 10;
				let chunkSize = 0;
				let chunk = 0;
				while (!request.resultBegin && chunk != 'RESULT_BEGIN') {
					chunkSize = parseInt(request.responseText.slice(startIndex, endIndex));
					startIndex 	= endIndex;
					endIndex 	= startIndex + chunkSize;
					chunk = request.responseText.slice(startIndex, endIndex);
					startIndex 	= endIndex;
					endIndex 	= startIndex + 10;
					request.position = startIndex;
				}

				responseText 	= window.atob(request.responseText.slice(request.position + 10).replace(/\s/g, ''));
				buffer 			= str2ab(responseText);
/*
	            bufferInflated 	= UZIP.inflate(buffer);
				console.log('bufferInflated.length : ' + bufferInflated.length);
*/
	            
/*
				returnBuffer 		= new Uint8Array(buffer);
				returnBufferSize 	= returnBuffer.byteLength;
*/
			
				//Write result buffer to file for debugging. Filename and MIME type are optional.
				//writetofile(returnBuffer, 'resultBuffer', 'application/octet-stream');
// 				md.results = parseresultsfrombuffer(md, bufferInflated, bufferInflated.length);
				md.results = parseresultsfrombuffer(md, buffer, buffer.length);
				
				// Let front end script handle changes to callout
				if (hasCallout) {
					callout.set('Success!', '<p>Request successful</p>');
				} else {
					solveButton.text(solveButtonText).prop('disabled', false);
				}
				successCallback();
			} catch (e) {
				console.log(e);
				if (vesl.strings.startsWith(responseText, 'Error')) {
					if (hasCallout) {
						callout.set(vesl.ERROR_HEADER_GENERAL, '<p>Something went wrong. ' + vesl.ERROR_CONTENT_TRY_AGAIN + '</p>');
					} else {
						solveButton.text('ISSM error!').prop('disabled', false);
					}
				} else {
					if (hasCallout) {
						callout.set(vesl.ERROR_HEADER_GENERAL, '<p>Something went wrong. ' + vesl.ERROR_CONTENT_TRY_AGAIN + '</p>');
					} else {
						solveButton.text('JS error!').prop('disabled', false);
					}
				}
				errorCallback();
			}
		}; //}}}
		
		request.responseType = 'text';
		
		/* Construct request */
		let npbuffer = this.str2ab(md.cluster.np.toString());
		npbuffer = UZIP.deflate(npbuffer);
		let nplength = new Uint32Array(1);
		nplength[0] = npbuffer.byteLength;
		
		let codeversionbuffer = this.str2ab(md.cluster.codeversion.toString());
		codeversionbuffer = UZIP.deflate(codeversionbuffer);
		let codeversionlength = new Uint32Array(1);
		codeversionlength[0] = codeversionbuffer.byteLength;
		
		let runtimenamebuffer = this.str2ab(runtimename);
		runtimenamebuffer = UZIP.deflate(runtimenamebuffer);
		let runtimenamelength = new Uint32Array(1);
		runtimenamelength[0] = runtimenamebuffer.byteLength;
		
		let namebuffer = this.str2ab(name);
		namebuffer = UZIP.deflate(namebuffer);
		let namelength = new Uint32Array(1);
		namelength[0] = namebuffer.byteLength;
		
		let toolkitsbuffer = this.str2ab(toolkitsstring);
		toolkitsbuffer = UZIP.deflate(toolkitsbuffer);
		let toolkitslength = new Uint32Array(1);
		toolkitslength[0] = toolkitsbuffer.byteLength;
		
		let solutionbuffer = this.str2ab(solutionstring);
		solutionbuffer = UZIP.deflate(solutionbuffer);
		let solutionlength = new Uint32Array(1);
		solutionlength[0] = solutionbuffer.byteLength;
		
		let binbuffer = new Uint8Array(fid.rawbuffer()); //seems that 16 bits length could be incompatible.
		binbuffer = UZIP.deflate(binbuffer);
		let binlength = new Uint32Array(1);
		binlength[0] = binbuffer.byteLength;
		
		let data = new Blob(
			[
				nplength,
				npbuffer,
				codeversionlength,
				codeversionbuffer,
				runtimenamelength,
				runtimenamebuffer,
				namelength,
				namebuffer,
				toolkitslength,
				toolkitsbuffer,
				solutionlength,
				solutionbuffer,
				binlength,
				binbuffer
			]
		);

		// Send request
		request.send(data);
		
		if (hasCallout) {
			callout.set('Connecting...', '');
		} else {
			solveButton.text('Connecting...').prop('disabled', true);
		}
	}; //}}}
	
	this.ab2str = function(buf) { //{{{
		return String.fromCharCode.apply(null, new Uint16Array(buf));
	}; //}}}
	
	this.str2ab = function(str) { //{{{
		let buf = new Uint8Array(str.length);
		
		for (let i = 0, strLen = str.length; i < strLen; i++) {
			buf[i] = str.charCodeAt(i);
		}
		
		return buf;
	}; //}}}
} 
