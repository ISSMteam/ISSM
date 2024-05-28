//BASAL FORCINGS class definition
//
//   Usage:
//      basalforcings=basalforcings();

function basalforcings(){
	//methods
	this.setdefaultparameters = function() {//{{{

	} // }}}
	this.disp = function(){ // {{{
		console.log(sprintf('   basal forcings parameters:'));

		fielddisplay(this,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]');
		fielddisplay(this,'floatingice_melting_rate','basal melting rate (positive if melting) [m/yr]');
		fielddisplay(this,'geothermalflux','geothermal heat flux [W/m^2]');

	} // }}}
    this.extrude = function(md) {//{{{
        this.groundedice_melting_rate=project3d(md,'vector',this.groundedice_melting_rate,'type','node','layer',1); 
        this.floatingice_melting_rate=project3d(md,'vector',this.floatingice_melting_rate,'type','node','layer',1); 
        this.geothermalflux=project3d(md,'vector',this.geothermalflux,'type','node','layer',1); //bedrock only gets geothermal flux
        return this;
    }//}}}
	this.classname = function(){ // {{{
		return "basalforcings";
	} // }}}
		this.initialize = function (md){ // {{{

			if (isNaN(this.groundedice_melting_rate)){
				this.groundedice_melting_rate=NewArrayFill(md.mesh.numberofvertices,0);
				console.log('      no basalforcings.groundedice_melting_rate specified: values set as zero');
			}

			if (isNaN(this.floatingice_melting_rate)){
				this.floatingice_melting_rate=NewArrayFill(md.mesh.numberofvertices,0);
				console.log('      no basalforcings.floatingice_melting_rate specified: values set as zero');
			}

		} // }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			if(ArrayAnyEqual(ArrayIsMember('MasstransportAnalysis',analyses),1)){
				if (!(solution=='TransientSolution' & md.trans.ismasstransport==0)){
					checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
					checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
				}
			}

			if(ArrayAnyEqual(ArrayIsMember('BalancethicknessAnalysis',analyses),1)){
				checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
			}
			if(ArrayAnyEqual(ArrayIsMember('ThermalAnalysis',analyses),1)){
				if (!(solution=='TransientSolution' & md.trans.isthermal==0)){
					checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
					checkfield(md,'fieldname','basalforcings.floatingice_melting_rate','NaN',1,'Inf',1,'timeseries',1);
					checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'timeseries',1,'>=',0);
				}
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			var yts=md.constants.yts;

			WriteData(fid,prefix,'name','md.basalforcings.model','data',1,'format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',this,'fieldname','floatingice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
			WriteData(fid,prefix,'object',this,'fieldname','geothermalflux','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties
	//{{{
	this.groundedice_melting_rate  = NaN;
	this.floatingice_melting_rate  = NaN;
	this.geothermalflux            = NaN;
	this.setdefaultparameters();
	//}}}
}
