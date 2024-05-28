//MASK class definition
//
//   Usage:
//      mask= new mask();

function mask () {
	//properties 
	// {{{
		this.ocean_levelset	= NaN;
		this.ice_levelset	= NaN;
		//}}}
	//methods 
		this.setdefaultparameters = function (){ //{{{
		} // }}}
		this.disp = function () { //{{{
			console.log(sprintf("   masks:")); 

			fielddisplay(this,'ocean_levelset','presence of ocean if < 0, coastline/grounding line if = 0, no ocean if > 0');
			fielddisplay(this,'ice_levelset','presence of ice if < 0, icefront position if = 0, no ice if > 0');
		} //}}}
		this.extrude = function(md) {//{{{
			this.ocean_levelset=project3d(md,'vector',this.ocean_levelset,'type','node');
			this.ice_levelset=project3d(md,'vector',this.ice_levelset,'type','node');
			return this;
		}//}}}
		this.classname = function () { //{{{
			return "mask";
		} //}}}
		this.checkconsistency = function(md,solution,analyses){ //{{{
			if (solution=='LoveSolution') return;

			checkfield(md,'fieldname','mask.ocean_levelset','size',[md.mesh.numberofvertices, 1]);
			checkfield(md,'fieldname','mask.ice_levelset'        ,'size',[md.mesh.numberofvertices, 1]);
			var isice=NewArrayFill(md.mesh.numberofvertices,0); 
			for(var i=0;i<md.mesh.numberofvertices;i++)if(md.mask.ice_levelset[i]<=0)isice[i]=1;
			if (ArraySum(isice)==0){
				console.log('no ice present in the domain');
			}
			if (ArrayMax(md.mask.ice_levelset)<0){
				console.log('no ice front provided');
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
            WriteData(fid,prefix,'object',this,'fieldname','ocean_levelset','format','DoubleMat','mattype',1);
            WriteData(fid,prefix,'object',this,'fieldname','ice_levelset','format','DoubleMat','mattype',1);
		}//}}}
		this.fix=function() { //{{{
		}//}}}
}
