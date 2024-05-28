async function solve(md, solutionstring, ...varargin) { //{{{
/*
SOLVE - apply solution sequence for this model

Usage:
	solve(md, solutionstring[, ...]);

where varargin is a list of paired arguments of string OR enums

Solution types available comprise:
- 'Stressbalance'			or 'sb'
- 'Masstransport'			or 'mt'
- 'Oceantransport'			or 'oceant'
- 'Thermal'					or 'th'
- 'Steadystate'				or 'ss'
- 'Transient'				or 'tr'
- 'Balancethickness'		or 'mc'
- 'BalancethicknessSoft'	or 'mcsoft'
- 'Balancevelocity'			or 'bv'
- 'BedSlope'				or 'bsl'
- 'SurfaceSlope'			or 'ssl'
- 'Hydrology'				or 'hy'
- 'DamageEvolution'			or 'da'
- 'Gia'						or 'gia'
- 'Love'					or 'lv'
- 'Esa'						or 'esa'
- 'Sampling'				or 'smp'
- 'Gmsh'

Extra options (these all need to be passed in via the third parameter, which is 
a rest parameter):
- loadonly    		: do not solve, only load results
- runtimename 		: true or false (default is true); makes name unique
- checkconsistency 	: true or false (default is true); checks consistency of model
- restart			: directory name (relative to the execution directory) where the restart file is located

Examples:
	md = solve(md, 'Stressbalance');
	md = solve(md, 'sb');
	
NOTE:
- We do not strictly need to return md as objects are passed by reference in 
JavaScript, but we do so to mirror MATLAB and Python APIs.

TODO:
- Refactor UI reporting structure so we do not have to check if it is defined
*/	
 
/*
	// Check that md exists and that it is a model
 	if (md === null || md === undefined || md.constructor.name !== 'model') {
 		throw new Error('md needs to be an instance of the model class');
 	}
*/
	
	if (typeof(solutionstring) !== 'string') {
		throw new Error('ISSM\'s solve function only accepts strings for solution sequences. Type help solve to get a list of supported solutions.');
	}

	// Recover and process solve options
	if (vesl.strings.strcmpi(solutionstring, 'sb') || vesl.strings.strcmpi(solutionstring, 'Stressbalance')) {
		solutionstring = 'StressbalanceSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'mt') || vesl.strings.strcmpi(solutionstring, 'Masstransport')) {
		solutionstring = 'MasstransportSolution';	
	} else if (vesl.strings.strcmpi(solutionstring, 'oceant') || vesl.strings.strcmpi(solutionstring, 'Oceantransport')) {
		solutionstring = 'OceantransportSolution';	
	} else if (vesl.strings.strcmpi(solutionstring, 'th') || vesl.strings.strcmpi(solutionstring, 'Thermal')) {
		solutionstring = 'ThermalSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'st') || vesl.strings.strcmpi(solutionstring, 'Steadystate')) {
		solutionstring = 'SteadystateSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'tr') || vesl.strings.strcmpi(solutionstring, 'Transient')) {
		solutionstring = 'TransientSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'mc') || vesl.strings.strcmpi(solutionstring, 'Balancethickness')) {
		solutionstring = 'BalancethicknessSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'mcsoft') || vesl.strings.strcmpi(solutionstring, 'BalancethicknessSoft')) {
		solutionstring = 'BalancethicknessSoftSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'bv') || vesl.strings.strcmpi(solutionstring, 'Balancevelocity')) {
		solutionstring = 'BalancevelocitySolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'bsl') || vesl.strings.strcmpi(solutionstring, 'BedSlope')) {
		solutionstring = 'BedSlopeSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'ssl') || vesl.strings.strcmpi(solutionstring, 'SurfaceSlope')) {
		solutionstring = 'SurfaceSlopeSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'hy') || vesl.strings.strcmpi(solutionstring, 'Hydrology')) {
		solutionstring = 'HydrologySolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'da') || vesl.strings.strcmpi(solutionstring, 'DamageEvolution')) {
		solutionstring = 'DamageEvolutionSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'gia') || vesl.strings.strcmpi(solutionstring, 'Gia')) {
		solutionstring = 'GiaSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'lv') || vesl.strings.strcmpi(solutionstring, 'Love')) {
		solutionstring = 'LoveSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'Esa')) {
		solutionstring = 'EsaSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'smp') || vesl.strings.strcmpi(solutionstring, 'Sampling')) {
		solutionstring = 'SamplingSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'gmsh')) {
		solutionstring = 'GmshSolution';
	} else if (vesl.strings.strcmpi(solutionstring, 'gmt')) {
		solutionstring = 'GmtSolution';
	} else {
		throw new Error('solutionstring ' + solutionstring + ' not supported!');
	}
	let options = new pairoptions(varargin, 'solutionstring', solutionstring);
	
	// Recover some fields
	md.priv.solution 	= solutionstring;
	let cluster 		= md.cluster;
	
	// NOTE: Batch scripts are not currently implemented
	let batch = 0; 
	if (options.getfieldvalue('batch', 'no') === 'yes') {
		batch = 1;
	}
	
	// Check model consistency
	if (options.getfieldvalue('checkconsistency', 'yes') === 'yes') {
		if (md.verbose.solution) {
			console.log('checking model consistency');
		}
		
		ismodelselfconsistent(md);
	}
	
	// If we are restarting, actually use the provided runtime name:
	restart = options.getfieldvalue('restart', '');

	// First, build a runtime name that is unique
	if (restart === 1) {
		// Leave the runtimename as is
	} else {
		if (restart !== '') {
			md.priv.runtimename = restart;
		} else {
			if (options.getfieldvalue('runtimename', true)) {
				let c = new Date().getTime();
				md.priv.runtimename = sprintf('%s-%g', md.miscellaneous.name, c);
			} else {
				md.priv.runtimename = md.miscellaneous.name;
			}
		}
	}

	// If running QMU analysis, some preprocessing of Dakota files using model fields needs to be carried out
	if (md.qmu.isdakota) {
		throw new Error("QMU not supported yet!");
		//md = preqmu(md, options);
	}
	
	// Do we load results only?
	if (options.getfieldvalue('loadonly', false)){
		loadresultsfromcluster(md);
		return;
	}

	/*
	Write all input arrays (as opposed to, under MATLAB/Python, input binary 
	files)
	
	NOTE: The JavaScript implementation diverges significantly from the 
		  MATLAB/Python APIs here.
	*/
	let fid = null; // bin file equivalent
	//TODO: FIND A BETTER WAY TO DO THIS! (IE, SYNC UP WRITEDATA AND HAVE A FULL DEMARSHALL/READMODEL IN PYTHON
	if (solutionstring === 'GmshSolution') {
		//open file for binary writing
		fid = new fileptr('mode','w');
	} else if (solutionstring === 'GmtSolution') {
		//open file for binary writing
		fid = new fileptr('mode','w');
		let prefix='md.mesh';
		WriteData(fid,prefix,'object',md.mesh,'fieldname','lat','format','DoubleMat','mattype',1);
		WriteData(fid,prefix,'object',md.mesh,'fieldname','long','format','DoubleMat','mattype',1);
	} else {
		// Marshall into a binary array (fid) all the fields of model
		fid = marshall(md); // bin file
	}
	let toolkitsstring = md.toolkits.ToolkitsFile(md.miscellaneous.name + '.toolkits'); // toolkits file equivalent

	if (cluster.classname() === 'local') {//{{{

		// We are running locally on the machine, using the ISSM module
		console.log('running issm locally');
		
		// Call ISSM
		let outputs = issm(fid, toolkitsstring, solutionstring, md.miscellaneous.name); 
		
		// Recover output
		let outputbuffer 		= outputs[0]; 
		let outputbuffersize 	= outputs[1];
			
		// Load results 
		md = loadresultsfrombuffer(md, outputbuffer, outputbuffersize); // TODO: Pass reporting construct to loadresultsfrombuffer
		
		// Call success callback
		if (vesl.helpers.isFunction(vesl.ui.reporting.success_callback)) {
			vesl.ui.reporting.success_callback();
		}
	//}}}
	} else { //{{{
		// We are running somewhere else on a computational server. Send the buffer to that server and retrieve output.
		console.log('running issm remotely');
		
		await cluster.uploadandrun(
			md, 
			fid, 
			toolkitsstring, 
			solutionstring, 
			md.miscellaneous.name, 
			md.priv.runtimename,
			options
		);/*
.catch(function(e) {
			if (vesl.helpers.isDefined(vesl.ui) && vesl.helpers.isDefined(vesl.ui.reporting) && vesl.helpers.isFunction(vesl.ui.reporting.error_callback)) {
				vesl.ui.reporting.error_callback(e);
			}
		}).catch(function(e) {
			// Handle unexpected errors (source: http://thecodebarbarian.com/async-await-error-handling-in-javascript.html)
			console.log(e);
		});
			
		if (vesl.helpers.isDefined(vesl.ui) && vesl.helpers.isDefined(vesl.ui.reporting) && vesl.helpers.isFunction(vesl.ui.reporting.success_callback)) {
			vesl.ui.reporting.success_callback(md);
		}
*/
		
		// Why is md undefined at vesl.ui.reporting.success_callback(md)? See issm-refactor
		
		return md;
	} //}}}
} //}}}
