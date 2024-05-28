//CALVING class definition
//
//   Usage:
//      calving=new calving();

function calving (){
	//methods
	this.setdefaultparameters = function(){// {{{

	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   Calving parameters:'));
		fielddisplay(this,'calvingrate','calving rate at given location [m/a]');

	}// }}}
    this.extrude = function(md) {//{{{
        this.calvingrate=project3d(md,'vector',this.calvingrate,'type','node');
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "calving";
	}// }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{
		//Early return
		if (solution!='TransientSolution' | md.trans.ismovingfront==0) return;

		checkfield(md,'fieldname','calving.calvingrate(1:md.mesh.numberofvertices,:)','>=',0,'timeseries',1,'NaN',1,'Inf',1);
	} //}}}
		this.marshall=function(md,prefix,fid) { //{{{
			var yts=md.constants.yts;
			WriteData(fid,prefix,'name','md.calving.law','data',1,'format','Integer');
			WriteData(fid,prefix,'object',this,'fieldname','calvingrate','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts,'scale',1./yts);
		}//}}}
		this.fix=function() { //{{{
			this.calvingrate=NullFix(this.calvingrate,NaN);
			this.meltingrate=NullFix(this.meltingrate,NaN);
		}//}}}
	//properties 
	// {{{

	this.calvingrate   = NaN;

	this.setdefaultparameters();
	//}}}
}
