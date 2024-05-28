//GEOMETRY class definition
//
//   Usage:
//      geometry=geometry();

function geometry(){
	//methods 
		this.setdefaultparameters = function (){ //{{{
		}// }}}
		this.disp = function () { //{{{
			console.log(sprintf("   Geometry parameters:"));

			fielddisplay(this,'surface','ice upper surface elevation [m]');
			fielddisplay(this,'thickness','ice thickness [m]');
			fielddisplay(this,'base','ice base elevation [m]');
			fielddisplay(this,'bed','bed elevation [m]');
		} //}}}
        this.extrude = function(md) {//{{{
            this.surface=project3d(md,'vector',this.surface,'type','node');
            this.thickness=project3d(md,'vector',this.thickness,'type','node');
            this.hydrostatic_ratio=project3d(md,'vector',this.hydrostatic_ratio,'type','node');
            this.base=project3d(md,'vector',this.base,'type','node');
            this.bed=project3d(md,'vector',this.bed,'type','node');
            return this;
        }//}}}
		this.classname = function () { //{{{
			return 'geometry';
		} //}}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			if ((solution=='TransientSolution' & md.trans.isgia) | (solution=='GiaSolution')){
				checkfield(md,'fieldname','geometry.thickness','timeseries',1,'NaN',1,'Inf',1,'>=',0);
			}
			else{
				checkfield(md,'fieldname','geometry.surface'  ,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','geometry.base'      ,'NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				checkfield(md,'fieldname','geometry.thickness','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1],'>',0);
				for(var i=0;i<md.mesh.numberofvertices;i++){
					if (Math.abs(md.geometry.thickness[i]-md.geometry.surface[i]+md.geometry.base[i])>Math.pow(10,9)){
						checkmessage(md,'equality thickness=surface-base violated');
						break;
					}
				}
				if (solution=='TransientSolution' & md.trans.isgroundingline){
					checkfield(md,'fieldname','geometry.bed','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices, 1]);
				}
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			let length_thickness=size(this.thickness,0);
			if (length_thickness==md.mesh.numberofvertices || length_thickness==md.mesh.numberofvertices+1) {
				WriteData(fid,prefix,'object',this,'fieldname','thickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts);
			} else if (length_thickness==md.mesh.numberofelements || length_thickness==md.mesh.numberofelements+1) {
				WriteData(fid,prefix,'object',this,'fieldname','thickness','format','DoubleMat','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts);
			} else {
				error('geometry thickness time series should be a vertex or element time series');
			}
			WriteData(fid,prefix,'object',this,'fieldname','surface','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','base','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','bed','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','hydrostatic_ratio','format','DoubleMat','mattype',1);
		}//}}}
		this.fix=function() { //{{{
			this.hydrostatic_ratio=NullFix(this.hydrostatic_ratio,NaN);
		}//}}}
	//properties 
	// {{{
		this.surface           = NaN;
		this.thickness         = NaN;
		this.base              = NaN;
		this.bed               = NaN;
		this.hydrostatic_ratio = NaN;
		this.setdefaultparameters();
		//}}}
}
