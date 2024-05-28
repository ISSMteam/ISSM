//FRONTAL FORCINGS Class definition
//
//   Usage:
//      frontalforcings=frontalforcings();

function frontalforcings(){
	//methods
    this.classname = function(){ // {{{
        return "frontalforcings";
    } // }}}
    this.extrude = function(md) {//{{{
		this.meltingrate=project3d(md,'vector',this.meltingrate,'type','node');
        return this;
    }//}}}
	this.setdefaultparameters = function(){// {{{
		this.meltingrate   = NaN;
	} // }}}
    this.checkconsistency = function(md,solution,analyses) { //{{{
		//Early return
		if (!solution=='TransientSolution' || md.transient.ismovingfront == 0) return;

		md = checkfield(md,'fieldname','frontalforcings.meltingrate','NaN',1,'Inf',1,'timeseries',1,'>=',0);

    } // }}}
	this.disp = function(){ // {{{
		console.log(sprintf('   Frontalforcings parameters:'));
		fielddisplay(this,'meltingrate','melting rate at given location [m/a]');
	} // }}}
    this.marshall=function(md,prefix,fid) { //{{{
		let yts=md.constants.yts;
		WriteData(fid,prefix,'name','md.frontalforcings.parameterization','data',1,'format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','meltingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
    }//}}}
    this.fix=function() { //{{{
    }//}}}
	//properties 
    // {{{
	this.meltingrate   = NaN;
	this.setdefaultparameters();
    // }}}
}
