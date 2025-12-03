//GROUNDINGLINE class definition
//
//   Usage:
//      groundingline=new groundingline();

function groundingline (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//Type of migration
		this.migration				= 'SubelementMigration';
		this.friction_interpolation	= 'SubelementFriction1';
		this.melt_interpolation		= 'NoMeltOnPartiallyFloating';
		this.intrusion_distance 	= 0;

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   grounding line migration parameters:'));
		fielddisplay(this,'migration',"type of grounding line migration: 'SubelementMigration','SoftMigration','AggressiveMigration','Contact' or 'None'");
		fielddisplay(this,'friction_interpolation',"type of friction interpolation forpartially floating elements: 'NoFrictionOnPartiallyFloating','SubelementFriction1' or 'SubelementFriction2'");
		fielddisplay(this,'melt_interpolation',"type of melt interpolation forpartially floating elements: 'NoMeltOnPartiallyFloating','FullMeltOnPartiallyFloating','SubelementMelt1','SubelementMelt2' or 'IntrusionMelt'");
		fielddisplay(this,'intrusion_distance','distance of seawater intrusion from grounding line [m]');

	}// }}}
	this.classname= function(){// {{{
		return "groundingline";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) {// {{{

			checkfield(md,'fieldname','groundingline.migration','values',['None', 'AggressiveMigration', 'SoftMigration', 'Contact', 'GroundingOnly']);
			checkfield(md,'fieldname','groundingline.friction_interpolation','values',['NoFrictionOnPartiallyFloating', 'SubelementFriction1', 'SubelementFriction2']);
			checkfield(md,'fieldname','groundingline.melt_interpolation','values',['NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating', 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt']);
			checkfield(md,'fieldname','groundingline.intrusion_distance','NaN',1,'Inf',1,'>=',0);

			if (this.migration !='None'){
				if (isNaN(md.geometry.bed)){
					md.checkmessage('requesting grounding line migration, but bathymetry is absent!');
				}
				for (var i=0;i<md.mesh.numberofvertices;i++){
					if(md.mask.groundedice_levelset[i]>0){
						md.checkmessage('base not equal to bed on grounded ice!');
						break;
					}
					if(md.geometry.bed[i] - md.geometry.base[i] > Math.pow(10,-9)){
						checkmessage(md,'bed superior to base on floating ice!');
						break;
					}
				}
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'data',this.migration,'name','md.groundingline.migration','format','String');
			WriteData(fid,prefix,'data',this.friction_interpolation,'name','md.groundingline.friction_interpolation','format','String');
			WriteData(fid,prefix,'data',this.melt_interpolation,'name','md.groundingline.melt_interpolation','format','String');
			WriteData(fid,prefix,'object',self,'fieldname','intrusion_distance','format','DoubleMat','mattype',1);

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.migration              = '';
	this.friction_interpolation = '';
	this.melt_interpolation     = '';
	this.intrusion_distance		= 0;
	this.setdefaultparameters();
	//}}}
}
