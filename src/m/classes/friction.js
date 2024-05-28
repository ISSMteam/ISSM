//FRICTION class definition
//
//   Usage:
//      friction=friction();

function friction (){
	//methods
	this.setdefaultparameters = function(){ // {{{

	} // }}}
	this.disp= function (){// {{{
		console.log(sprintf('Basal shear stress parameters: Sigma_b = coefficient^2 * Neff ^r * |u_b|^(s-1) * u_b\n(effective stress Neff=rho_ice*g*thickness+rho_water*g*bed, r=q/p and s=1/p)'));
		fielddisplay(this,'coefficient','friction coefficient [SI]');
		fielddisplay(this,'p','p exponent');
		fielddisplay(this,'q','q exponent');
		fielddisplay(this,'effective_pressure','Effective Pressure for the forcing if not coupled [Pa]');
		fielddisplay(this,'coupling','Coupling flag: 0 for default, 1 for forcing(provide md.friction.effective_pressure)  and 2 for coupled(not implemented yet)');
		fielddisplay(this,'effective_pressure_limit','Neff do not allow to fall below a certain limit: effective_pressure_limit*rho_ice*g*thickness (default 0)');
	} // }}}
	this.extrude = function(md) {//{{{
		this.coefficient = project3d(md, 'vector', this.coefficient, 'type', 'node', 'layer', 1);
		this.p = project3d(md, 'vector', this.p, 'type', 'element');
		this.q = project3d(md, 'vector', this.q, 'type', 'element');
		switch (this.coupling) {
			case 0:
			case 1:
				this.effective_pressure=project3d(md, 'vector', this.effective_pressure, 'type', 'node', 'layer', 1);
				break;
			case 2:
				console.error('not implemented yet');
				break;
			default:	
				console.error('not supported yet');
		}
		return this;
	}//}}}
	this.classname= function (){// {{{
		return "friction";
	} // }}}
		this.checkconsistency = function(md,solution,analyses){ //{{{

			//Early return
			if ((!ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1)) & (!ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1))){
				return; 
			}
			checkfield(md,'fieldname','friction.coefficient','timeseries',1,'NaN',1,'Inf',1);
			checkfield(md,'fieldname','friction.q','NaN',1,'Inf',1,'size',[md.mesh.numberofelements ,1]);
			checkfield(md,'fieldname','friction.p','NaN',1,'Inf',1,'size',[md.mesh.numberofelements ,1]);
			checkfield(md,'fieldname','friction.coupling','numel',[1],'values',[0, 1, 2]);
			checkfield(md,'fieldname','friction.effective_pressure_limit','numel',[1],'>=',0);
			switch (this.coupling) {
				case 0:
				case 1:
					checkfield(md,'fieldname','friction.effective_pressure','NaN',1,'Inf',1,'timeseries',1);
					break;
				case 2:
					console.error('not implemented yet');
					break;
				default:	
					console.error('not supported yet');
			}

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			var yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.friction.law','data',1,'format','Integer');
			let mattype,tsl;
			if ((size(this.coefficient,0)==md.mesh.numberofvertices | size(this.coefficient,1)==md.mesh.numberofvertices+1)) {
				mattype=1;
				tsl = md.mesh.numberofvertices;
			} else {
				mattype=2;
				tsl = md.mesh.numberofelements;
			}
			WriteData(fid,prefix,'object',this,'fieldname','coefficient','format','DoubleMat','mattype',mattype,'timeserieslength',tsl+1,'yts',md.constants.yts);
			WriteData(fid,prefix,'object',this,'fieldname','p','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'object',this,'fieldname','q','format','DoubleMat','mattype',2);
			WriteData(fid,prefix,'class','friction','object',this,'fieldname','coupling','format','Integer');
			WriteData(fid,prefix,'object',this,'class','friction','fieldname','effective_pressure_limit','format','Double');
			switch (this.coupling) {
				case 0:
					break;
				case 1:
					WriteData(fid,prefix,'class','friction','object',this,'fieldname','effective_pressure','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
					break;
				case 2:
					console.error('not implemented yet');
					break;
				default:
					console.error('not supported yet');		
			}
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	//{{{
	this.coefficient			  = NaN;
	this.p						  = NaN;
	this.q						  = NaN;
	this.coupling				  = 0;
	this.effective_pressure 	  = NaN;
	this.effective_pressure_limit = 0;
	this.setdefaultparameters();
	//}}}
}
