//LEVELSET class definition
//
//   Usage:
//      levelset=new levelset();

function levelset (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//stabilization = 1 by default
		this.stabilization		= 1;
		this.reinit_frequency	= 5;
		this.kill_icebergs      = 1;
		this.migration_max      = 1e12; //No need for general cases, unless specified

		//Linear elements by default
		this.fe='P1';
	
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   Level-set parameters:'));
		fielddisplay(this,'stabilization','0: no, 1: artificial_diffusivity, 2: streamline upwinding');
		fielddisplay(this,'spclevelset','Levelset constraints (NaN means no constraint)');
		fielddisplay(this,'reinit_frequency','Amount of time steps after which the levelset function in re-initialized (NaN: no re-initialization).');
		fielddisplay(this,'kill_icebergs','remove floating icebergs to prevent rigid body motions (1: true, 0: false)');
		fielddisplay(this,'migration_max','maximum allowed migration rate (m/a)');
		fielddisplay(this,'fe','Finite Element type: "P1" (default), or "P2"');

	}// }}}
    this.extrude = function(md) {//{{{
        this.spclevelset=project3d(md,'vector',this.spclevelset,'type','node');
        return this;
    }//}}}
	this.classname= function(){// {{{
		return "levelset";
	}// }}}
	this.checkconsistency = function(md,solution,analyses) { // {{{
		//Early return
		if (solution!='TransientSolution' | md.trans.ismovingfront==0) return;

		checkfield(md,'fieldname','levelset.spclevelset','Inf',1,'timeseries',1);
		checkfield(md,'fieldname','levelset.stabilization','values',[0,1,2,5]);
		checkfield(md,'fieldname','levelset.kill_icebergs','numel',1,'values',[0, 1]);
		checkfield(md,'fieldname','levelset.migration_max','numel',1,'NaN',1,'Inf',1,'>',0);
		checkfield(md,'fieldname','levelset.fe','values',['P1','P2']);
	} //}}}
	this.marshall=function(md,prefix,fid) { //{{{

		let yts=md.constants.yts;

		WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','spclevelset','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
		WriteData(fid,prefix,'object',this,'fieldname','reinit_frequency','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','kill_icebergs','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','migration_max','format','Double','scale',1/yts);
		WriteData(fid,prefix,'object',this,'fieldname','fe','format','String');
	}//}}}
		this.fix=function() { //{{{
			this.spclevelset=NullFix(this.spclevelset,NaN);
		}//}}}
	//properties 
	// {{{

	this.stabilization		= 0;
	this.spclevelset		= NaN;
	this.reinit_frequency	= NaN;
	this.kill_icebergs     	= 0;
	this.migration_max      = 0.;
	this.fe              	= 'P1';

	this.setdefaultparameters();
	//}}}
}
