//FOURIERLOVE class definition
//
//   Usage:
//      md.love=fourierlove();

function fourierlove (){
	//methods
	this.extrude = function(md) { // {{{
		return this;
	} // }}}
	this.setdefaultparameters = function() { // {{{
		//we setup an elastic love number computation by default.
			this.nfreq=1; 
			this.frequencies=[0]; //Hz
			this.sh_nmax=256; // .35 degree, 40 km at the equator.
			this.sh_nmin=1;
			// work on matlab script for computing g0 for given Earth's structure. 
			this.g0=9.81; // m/s^2; 
			this.r0=6371*1e3; //m;
			this.mu0=1e11; // Pa
			this.Gravitational_Constant=6.67259e-11; // m^3 kg^-1 s^-2
			this.allow_layer_deletion=1;
			this.underflow_tol=1e-16; //threshold of deep to surface love number ratio to trigger the deletion of layer 
			this.integration_steps_per_layer=100;
			this.istemporal=0;
			this.n_temporal_iterations=8;
			this.time=[0]; //s
			this.love_kernels=0; 
			this.forcing_type = 11; // surface loading
			this.inner_core_boundary=1;
			this.core_mantle_boundary=2;
			this.complex_computation=0;
	} // }}}
	this.disp = function() { // {{{
		fielddisplay(this,'nfreq','number of frequencies sampled (default 1, elastic) [Hz]');
		fielddisplay(this,'frequencies','frequencies sampled (convention defaults to 0 for the elastic case) [Hz]');
		fielddisplay(this,'sh_nmax','maximum spherical harmonic degree (default 256, .35 deg, or 40 km at equator)');
		fielddisplay(this,'sh_nmin','minimum spherical harmonic degree (default 1)');
		fielddisplay(this,'g0','adimensioning constant for gravity (default 10) [m/s^2]');
		fielddisplay(this,'r0','adimensioning constant for radius (default 6371*10^3) [m]');
		fielddisplay(this,'mu0','adimensioning constant for stress (default 10^11) [Pa]');
		fielddisplay(this,'Gravitational_Constant','Newtonian constant of gravitation (default 6.67259e-11 [m^3 kg^-1 s^-2])');
		fielddisplay(this,'allow_layer_deletion','allow for migration of the integration boundary with increasing spherical harmonics degree (default 1)');
		fielddisplay(this,'underflow_tol','threshold of deep to surface love number ratio to trigger the deletion of layers (default 2.2204460492503131E-016)');
		fielddisplay(this,'integration_steps_per_layer','number of radial steps to propagate the yi system from the bottom to the top of each layer (default 100)');
		fielddisplay(this,'istemporal','1 for time-dependent love numbers, 0 for frequency-dependent or elastic love numbers (default 0)', 'If 1: use fourierlove function build_frequencies_from_time to meet consistency');
		fielddisplay(this,'n_temporal_iterations','max number of iterations in the inverse Laplace transform. Also the number of spectral samples per time step requested (default 8)');
		fielddisplay(this,'time','time vector for deformation if istemporal (default 0) [s]');
		fielddisplay(this,'love_kernels','compute love numbers at depth? (default 0)');
		fielddisplay(this,'forcing_type','integer indicating the nature and depth of the forcing for the Love number calculation (default 11) :','1:  Inner core boundary -- Volumic Potential','2:  Inner core boundary -- Pressure','3:  Inner core boundary -- Loading','4:  Inner core boundary -- Tangential traction','5:  Core mantle boundary -- Volumic Potential','6:  Core mantle boundary -- Pressure','7:  Core mantle boundary -- Loading','8:  Core mantle boundary -- Tangential traction','9:  Surface -- Volumic Potential','10: Surface -- Pressure','11: Surface -- Loading','12: Surface -- Tangential traction ');
		fielddisplay(this,'inner_core_boundary','interface index in materials.radius locating forcing. Only used for forcing_type 1--4 (default 1)');
		fielddisplay(this,'core_mantle_boundary','interface index in materials.radius locating forcing. Only used for forcing_type 5--8 (default 2)'); 

	} // }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{

		if (ArrayAnyEqual(ArrayIsMember('LoveAnalysis',analyses),1)) return; 

		checkfield(md,'fieldname','love.nfreq','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.frequencies','NaN',1,'Inf',1,'numel',md.love.nfreq);
		checkfield(md,'fieldname','love.sh_nmax','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.sh_nmin','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.g0','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.r0','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.mu0','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.Gravitational_Constant','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.allow_layer_deletion','values',[0, 1]);
		checkfield(md,'fieldname','love.underflow_tol','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.integration_steps_per_layer','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','love.love_kernels','values',[0, 1]);
		checkfield(md,'fieldname','love.forcing_type','NaN',1,'Inf',1,'numel',1,'>',0, '<=', 12);
		checkfield(md,'fieldname','love.complex_computation','NaN',1,'Inf',1,'numel',1,'values',[0, 1]);

		checkfield(md,'fieldname','love.istemporal','values',[0, 1]);

		if (md.love.istemporal==1){
			checkfield(md,'fieldname','love.n_temporal_iterations','NaN',1,'Inf',1,'numel',1,'>',0);
			checkfield(md,'fieldname','love.time','NaN',1,'Inf',1,'numel',md.love.nfreq/2/md.love.n_temporal_iterations);
		}
		if (md.love.sh_nmin<=1 && md.love.forcing_type==9 || md.love.forcing_type==5 || md.love.forcing_type==1) {
			throw 'Degree 1 not supported for Volumetric Potential forcing. Use sh_min>=2 for this kind of calculation.';
		}

		//need 'litho' material: 
		console.log('md.fourierlove check consistency only paritally implemented for litho material');
		/*
		if ~isa(md.materials,'materials') | ~sum(strcmpi(md.materials.nature,'litho'))
			error('Need a "litho" material to run a Fourier Love number analysis');
		end

		mat=find(strcmpi(md.materials.nature,'litho'));
		if (md.love.forcing_type<=4) {
			checkfield(md,'fieldname','love.inner_core_boundary','NaN',1,'Inf',1,'numel',1,'>',0, '<=', md.materials(mat).numlayers);
		} else if (md.love.forcing_type<=8) {
			checkfield(md,'fieldname','love.core_mantle_boundary','NaN',1,'Inf',1,'numel',1,'>',0, '<=', md.materials(mat).numlayers);
		} */
	} // }}}
	this.marshall = function(md,prefix,fid) { // {{{
	
		WriteData(fid,prefix,'object',this,'fieldname','nfreq','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','frequencies','format','DoubleMat','mattype',3);
		WriteData(fid,prefix,'object',this,'fieldname','sh_nmax','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','sh_nmin','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','g0','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','r0','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','mu0','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','Gravitational_Constant','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','allow_layer_deletion','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','underflow_tol','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','integration_steps_per_layer','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','istemporal','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','n_temporal_iterations','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','complex_computation','format','Boolean');
		//note: no need to marshall the time vector, we have frequencies
		WriteData(fid,prefix,'object',this,'fieldname','love_kernels','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','forcing_type','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','inner_core_boundary','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','core_mantle_boundary','format','Integer');

	} // }}}
	//properties 
	// {{{
	this.nfreq                		= NaN;
	this.frequencies          		= NaN;
	this.sh_nmax              		= NaN;
	this.sh_nmin              		= NaN;
	this.g0                   		= NaN; 
	this.r0                   		= NaN; 
	this.mu0                  		= NaN;
	this.Gravitational_Constant 	= 0;
	this.allow_layer_deletion 		= NaN;
	this.underflow_tol              = 0;
	this.integration_steps_per_layer= 0;
	this.istemporal		   			= 0;
	this.n_temporal_iterations	  	= 0;
	this.time			            = 0;
	this.love_kernels 				= NaN;
	this.forcing_type         		= NaN;
	this.inner_core_boundary	    = 0;
	this.core_mantle_boundary	    = 0;
	this.complex_computation        = 0;
	
	//set defaults
	this.setdefaultparameters();
	// }}}
}
