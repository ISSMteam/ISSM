//SLR class definition
//
//   Usage:
//      slr=slr();

function slr(){
	//methods
		this.setdefaultparameters = function (){ //{{{

		//Convergence criterion: absolute, relative and residual
		this.reltol=0.01; // 1 per cent
		this.abstol=NaN;  //default

		//maximum of non-linear iterations.
		this.maxiter=5;

		//computational flags:
		this.rigid=1;
		this.elastic=1;
		this.rotation=0;
		this.ocean_area_scaling=0;

		//tidal love numbers:
		this.tide_love_h=0.6149; //degree 2
		this.tide_love_k=0.3055; //degree 2

		//secular fluid love number:
		this.fluid_love=0.942;

		//moment of inertia:
		this.equatorial_moi=8.0077*10^37; // [kg m^2]
		this.polar_moi		 =8.0345*10^37; // [kg m^2]

		// mean rotational velocity of earth:
		this.angular_velocity=7.2921*10^-5; // [s^-1]

		//numerical discretization accuracy
		this.degacc=.01;

		//steric:
		this.steric_rate=0;


		//output default:
		this.requested_outputs=['default'];

		//transitions should be a cell array of vectors:
		this.transitions=[];

		}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			//Early return
			if (ArrayAnyEqual(ArrayIsMember('SealevelriseAnalysis',analyses),0) || ArrayAnyEqual(ArrayIsMember('TransientSolution',analyses),0) && !md.transient.isslr) {
				return;
			}

			md = checkfield(md,'fieldname','slr.deltathickness','NaN',1,'Inf',1,'size',[md.mesh.numberofelements, 1]);
			md = checkfield(md,'fieldname','slr.sealevel','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','slr.love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.love_k','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.love_l','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.tide_love_h','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.tide_love_k','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.fluid_love','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.equatorial_moi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.polar_moi','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.angular_velocity','NaN',1,'Inf',1);
			md = checkfield(md,'fieldname','slr.reltol','size',[1, 1]);
			md = checkfield(md,'fieldname','slr.abstol','size',[1, 1]);
			md = checkfield(md,'fieldname','slr.steric_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			md = checkfield(md,'fieldname','slr.maxiter','size',[1, 1],'>=',1);
			md = checkfield(md,'fieldname','slr.degacc','size',[1, 1],'>=',1e-10);
			md = checkfield(md,'fieldname','slr.requested_outputs','stringrow',1);

			//check that love numbers are provided at the same level of accuracy:
			if (this.love_h.length != this.love_k.length || this.love_h.length != this.love_l.length){
				throw Error('slr error message: love numbers should be provided at the same level of accuracy');
			}

		} // }}}
		this.defaultoutputs = function(md){ // {{{
			return ['Sealevel'];
		}//}}}
		this.classname= function(){// {{{
			return "slr";
		}// }}}
		this.disp= function(){// {{{

			console.log(sprintf('   Sealevelrise solution parameters:'));

		fielddisplay(this,'deltathickness','thickness change (main loading of the slr solution core [m]');
		fielddisplay(this,'sealevel','current sea level (prior to computation) [m]');
		fielddisplay(this,'reltol','sea level rise relative convergence criterion, (default, NaN: not applied)');
		fielddisplay(this,'abstol','sea level rise absolute convergence criterion, NaN: not applied');
		fielddisplay(this,'maxiter','maximum number of nonlinear iterations');
		fielddisplay(this,'love_h','load Love number for radial displacement');
		fielddisplay(this,'love_k','load Love number for gravitational potential perturbation');
		fielddisplay(this,'love_l','load Love number for horizontal displacements');
		fielddisplay(this,'tide_love_h','tidal love number (degree 2)');
		fielddisplay(this,'tide_love_k','tidal love number (degree 2)');
		fielddisplay(this,'fluid_love','secular fluid Love number');
		fielddisplay(this,'equatorial_moi','mean equatorial moment of inertia [kg m^2]');
		fielddisplay(this,'polar_moi','polar moment of inertia [kg m^2]');
		fielddisplay(this,'angular_velocity','mean rotational velocity of earth [per second]');
		fielddisplay(this,'rigid','rigid earth graviational potential perturbation');
		fielddisplay(this,'elastic','elastic earth graviational potential perturbation');
		fielddisplay(this,'rotation','rotational earth potential perturbation');
		fielddisplay(this,'ocean_area_scaling','correction for model representation of ocean area [default: No correction]');
		fielddisplay(this,'degacc',"accuracy (default .01 deg) for numerical discretization of the Green's functions");
		fielddisplay(this,'transitions','indices into parts of the mesh that will be icecaps');
		fielddisplay(this,'requested_outputs','additional outputs requested');
		fielddisplay(this,'steric_rate','rate of steric ocean expansion [mm/yr]');
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{
			console.log('WARNING: NOT MARHSALLING SLR FOR NOW.');
			return;
			WriteData(fid,prefix,'object',this,'fieldname','deltathickness','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'fieldname','sealevel','mattype',1,'format','DoubleMat','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','reltol','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','abstol','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','maxiter','format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','love_h','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','love_k','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','love_l','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','tide_love_h','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','tide_love_k','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','fluid_love','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','equatorial_moi','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','polar_moi','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','angular_velocity','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','rigid','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','elastic','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','rotation','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','ocean_area_scaling','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','steric_rate','format','DoubleMat','mattype',1,'scale',1e-3/md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','degacc','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','transitions','format','MatArray');

			//process requested outputs
			var outputs = this.requested_outputs;
			for (var i=0;i<outputs.length;i++){
				if (outputs[i] == 'default') {
					outputs.splice(i,1);
					var newoutputs=this.defaultoutputs(md);
					for (var j=0;j<newoutputs.length;j++) outputs.push(newoutputs[j]);
				}
			}
			WriteData(fid,prefix,'data',outputs,'name','md.slr.requested_outputs','format','StringArray');
		}//}}}
		this.fix=function() { //{{{
			this.deltathickness=NullFix(this.deltathickness,NaN);
			this.sealevel=NullFix(this.sealevel,NaN);
			this.maxiter=NullFix(this.maxiter,NaN);
			this.reltol=NullFix(this.reltol,NaN);
			this.abstol=NullFix(this.abstol,NaN);
			this.love_h=NullFix(this.love_h,NaN);
			this.love_k=NullFix(this.love_k,NaN);
			this.love_l=NullFix(this.love_l,NaN);
			this.tide_love_h=NullFix(this.tide_love_h,NaN);
			this.tide_love_k=NullFix(this.tide_love_k,NaN);
			this.fluid_love=NullFix(this.fluid_love,NaN);
			this.equatorial_moi=NullFix(this.equatorial_moi,NaN);
			this.polar_moi=NullFix(this.polar_moi,NaN);
			this.angular_velocity=NullFix(this.angular_velocity,NaN);
			this.rigid=NullFix(this.rigid,NaN);
			this.elastic=NullFix(this.elastic,NaN);
			this.rotation=NullFix(this.rotation,NaN);
			this.ocean_area_scaling=NullFix(this.ocean_area_scaling,NaN);
			this.steric_rate=NullFix(this.steric_rate,NaN);
			this.degacc=NullFix(this.degacc,NaN);
		}//}}}
	//properties
	//{{{
	this.deltathickness = NaN;
	this.sealevel       = NaN;
	this.maxiter        = 0;
	this.reltol         = 0;
	this.abstol         = 0;
	this.love_h         = 0; //provided by PREM model
	this.love_k         = 0; //idam
	this.love_l         = 0; //idam
	this.tide_love_h    = 0;
	this.tide_love_k    = 0;
	this.fluid_love	  = 0;
	this.equatorial_moi	= 0;
	this.polar_moi			= 0;
	this.angular_velocity = 0;
	this.rigid          = 0;
	this.elastic        = 0;
	this.rotation       = 0;
	this.ocean_area_scaling = 0;
	this.steric_rate    = 0; //rate of ocean expansion from steric effects.
	this.degacc         = 0;
	this.requested_outputs = [];
	this.transitions    = [];
	this.setdefaultparameters();
	//}}}
}
