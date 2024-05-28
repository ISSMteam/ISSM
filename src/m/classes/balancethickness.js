//BALANCETHICKNESS class definition
//
//   Usage:
//      balancethickness=new balancethickness();

function balancethickness (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//Type of stabilization used
		this.stabilization=1;

	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   balance thickness solution parameters:'));

		fielddisplay(this,'spcthickness','thickness constraints (NaN means no constraint) [m]');
		fielddisplay(this,'thickening_rate','ice thickening rate used in the mass conservation (dh/dt) [m/yr]');
		fielddisplay(this,'stabilization',"0: None, 1: SU, 2: SSA's artificial diffusivity, 3:DG");

	}// }}}
	this.classname= function(){// {{{
		return "balancethickness";

	}// }}}
		this.checkconsistency = function(md,solution,analyses){ // {{{
			//Early return
			if (solution!='BalancethicknessSolution')return;

			checkfield(md,'fieldname','balancethickness.spcthickness');
			checkfield(md,'fieldname','balancethickness.thickening_rate','size',[md.mesh.numberofvertices ,1],'NaN',1,'Inf',1);
			checkfield(md,'fieldname','balancethickness.stabilization','size',[1, 1],'values',[0, 1, 2 ,3]);
			//checkfield(md,'fieldname','balancethickness.omega','size',[md.mesh.numberofvertices ,1],'NaN',1,'Inf',1,'>=',0);
		} //}}}
		this.marshall=function(md,prefix,fid) { //{{{

			var yts=md.constants.yts;

			WriteData(fid,prefix,'object',this,'fieldname','spcthickness','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','thickening_rate','format','DoubleMat','mattype',1,'scale',1/yts);
			WriteData(fid,prefix,'object',this,'fieldname','stabilization','format','Integer');

			WriteData(fid,prefix,'object',this,'fieldname','slopex','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','slopey','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','omega','format','DoubleMat','mattype',1);

		}//}}}
		this.fix=function() { //{{{
			this.spcthickness=NullFix(this.spcthickness,NaN);
			this.thicknening_rate=NullFix(this.thicknening_rate,NaN);
			this.omega=NullFix(this.omega,NaN);
		}//}}}
	//properties 
	// {{{
	this.spcthickness      = NaN;
	this.thickening_rate   = NaN;
	this.stabilization     = 0;

	this.omega             = NaN;
	this.slopex            = NaN;
	this.slopey            = NaN;
	this.setdefaultparameters();
	//}}}
}
