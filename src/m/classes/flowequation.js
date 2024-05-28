//FLOWEQUATION class definition
//
//   Usage:
//      flowequation=new flowequation();

function flowequation (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//P1 for SSA
		this.fe_SSA= 'P1';

		//P1 for HO
		this.fe_HO= 'P1';

		//MINI condensed element for FS by default
		this.fe_FS = 'MINIcondensed';
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   flow equation parameters:'));

		fielddisplay(this,'isSIA','is the Shallow Ice Approximation (SIA) used ?');
		fielddisplay(this,'isSSA','is the Shelfy-Stream Approximation (SSA) used ?');
		fielddisplay(this,'isL1L2','is the L1L2 approximation used ?');
		fielddisplay(this,'isMOLHO','is the MOno-Layer Higher-Order (MOLHO) approximation used?');
		fielddisplay(this,'isHO','is the Higher-Order (HO) approximation used ?');
		fielddisplay(this,'isFS','are the Full-FS (FS) equations used ?');
		fielddisplay(this,'isNitscheBC','is weakly imposed condition used?');
		fielddisplay(this,'FSNitscheGamma','Gamma value for the Nitsche term (default: 1e6)');
		fielddisplay(this,'fe_SSA',"Finite Element for SSA  'P1', 'P1bubble' 'P1bubblecondensed' 'P2'");
		fielddisplay(this,'fe_HO', "Finite Element for HO   'P1' 'P1bubble' 'P1bubblecondensed' 'P1xP2' 'P2xP1' 'P2'");
		fielddisplay(this,'fe_FS', "Finite Element for FS   'P1P1' (debugging only) 'P1P1GLS' 'MINIcondensed' 'MINI' 'TaylorHood' 'XTaylorHood'");
		fielddisplay(this,'vertex_equation','flow equation for each vertex');
		fielddisplay(this,'element_equation','flow equation for each element');
		fielddisplay(this,'borderSSA',"vertices on SSA's border (for tiling)");
		fielddisplay(this,'borderHO',"vertices on HO's border (for tiling)");
		fielddisplay(this,'borderFS',"vertices on FS' border (for tiling)");

	}// }}}
	this.classname= function(){// {{{
		return "flowequation";

	}// }}}
	this.extrude = function(md) {//{{{
		this.element_equation=project3d(md,'vector',this.element_equation,'type','element');
		this.vertex_equation=project3d(md,'vector',this.vertex_equation,'type','node');
		this.borderSSA=project3d(md,'vector',this.borderSSA,'type','node');
		this.borderHO=project3d(md,'vector',this.borderHO,'type','node');
		this.borderFS=project3d(md,'vector',this.borderFS,'type','node');
		return this;
    }//}}}
		this.checkconsistency = function(md,solution,analyses) {//{{{

			//Early return
			if ( ((!ArrayAnyEqual(ArrayIsMember('StressbalanceAnalysis',analyses),1)) & (!ArrayAnyEqual(ArrayIsMember('StressbalanceSIAAnalysis',analyses),1))) | 
					(solution=='TransientSolution' & md.trans.isstressbalance==0)
			   ) return ;

			checkfield(md,'fieldname','flowequation.isSIA','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isSSA','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isL1L2','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isMOLHO','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isHO','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isFS','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.isNitscheBC','numel',[1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.FSNitscheGamma','numel',[1], '>=', 0);
			checkfield(md,'fieldname','flowequation.fe_SSA','values',['P1','P1bubble','P1bubblecondensed','P2','P2bubble']);
			checkfield(md,'fieldname','flowequation.fe_HO' ,'values',['P1','P1bubble','P1bubblecondensed','P1xP2','P2xP1','P2','P2bubble','P1xP3','P2xP4']);
			checkfield(md,'fieldname','flowequation.fe_FS' ,'values',['P1P1','P1P1GLS','MINIcondensed','MINI','TaylorHood','LATaylorHood','XTaylorHood','OneLayerP4z','CrouzeixRaviart','LACrouzeixRaviart']);
			checkfield(md,'fieldname','flowequation.augmented_lagrangian_r','numel',[1],'>=',0.);
			checkfield(md,'fieldname','flowequation.augmented_lagrangian_rlambda','numel',[1],'>=',0.);
			checkfield(md,'fieldname','flowequation.augmented_lagrangian_rhop','numel',[1],'>=',0.);
			checkfield(md,'fieldname','flowequation.augmented_lagrangian_rholambda','numel',[1],'>=',0.);
			checkfield(md,'fieldname','flowequation.XTH_theta','numel',[1],'>=',0.,'<',0.5);
			checkfield(md,'fieldname','flowequation.borderSSA','size',[md.mesh.numberofvertices, 1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.borderHO','size',[md.mesh.numberofvertices, 1],'values',[0, 1]);
			checkfield(md,'fieldname','flowequation.borderFS','size',[md.mesh.numberofvertices, 1],'values',[0, 1]);
			if (md.mesh.domaintype() == '2Dhorizontal'){
				checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices, 1],'values',[1,2]);
				checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements, 1],'values',[1,2]);
			}
			else if (md.mesh.domaintype() == '3Dsurface'){
				checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices, 1],'values',[1,2]);
				checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements, 1],'values',[1,2]);
			}
			else if (md.mesh.domaintype() =='2Dvertical'){
				checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices, 1],'values',[2,4,5]);
				checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements, 1],'values',[2,4,5]);
			}
			else if (md.mesh.domaintype() =='3D'){
				checkfield(md,'fieldname','flowequation.vertex_equation','size',[md.mesh.numberofvertices, 1],'values',[0,1,2,3,4,5,6,7,8]);
				checkfield(md,'fieldname','flowequation.element_equation','size',[md.mesh.numberofelements, 1],'values',[0,1,2,3,4,5,6,7,8]);
			}
			else throw Error('Case not supported yet');
			
			if (!(this.isSIA | this.isSSA | this.isL1L2 | this.isMOLHO | this.isHO | this.isFS)){
				checkmessage(md,['no element types set for this model']);
			}
			if(ArrayAnyEqual(ArrayIsMember('StressbalanceSIAAnalysis', analyses),1)){
				if (ArrayAnyEqual(this.element_equation,1)){
					if(this.vertex_equation & ArrayAnyBelowStrict(md.mask.groundedice_levelset)){
						console.log(sprintf("\n !!! Warning: SIA's model is not consistent on ice shelves !!!\n"));
					}
				}
			}
		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','isSIA','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isSSA','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isL1L2','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isMOLHO','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isHO','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isFS','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','isNitscheBC','format','Boolean');
			WriteData(fid,prefix,'object',this,'fieldname','FSNitscheGamma','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','fe_SSA','data',this.fe_SSA,'format','String');
			WriteData(fid,prefix,'object',this,'fieldname','fe_HO','data',this.fe_HO,'format','String');
			WriteData(fid,prefix,'object',this,'fieldname','fe_FS','data',this.fe_FS,'format','String');

			WriteData(fid,prefix,'object',this,'fieldname','augmented_lagrangian_r','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','augmented_lagrangian_rhop','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','augmented_lagrangian_rlambda','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','augmented_lagrangian_rholambda','format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','XTH_theta','data',this.XTH_theta ,'format','Double');
			WriteData(fid,prefix,'object',this,'fieldname','borderSSA','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','borderHO','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'object',this,'fieldname','borderFS','format','DoubleMat','mattype',1);

			//convert approximations to integers 
			WriteData(fid,prefix,'data',this.vertex_equation,'name','md.flowequation.vertex_equation','format','DoubleMat','mattype',1);
			WriteData(fid,prefix,'data',this.element_equation,'name','md.flowequation.element_equation','format','DoubleMat','mattype',2);

		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.isSIA                          = 0;
	this.isSSA                          = 0;
	this.isL1L2                         = 0;
	this.isMOLHO                         = 0;
	this.isHO                           = 0;
	this.isFS                           = 0;
	this.isNitscheBC                    = 0;
	this.FSNitscheGamma                 = 1e6;
	this.fe_SSA                         = '';
	this.fe_HO                          = '';
	this.fe_FS                          = '';
	this.augmented_lagrangian_r         = 1.;
	this.augmented_lagrangian_rhop      = 1.;
	this.augmented_lagrangian_rlambda   = 1.;
	this.augmented_lagrangian_rholambda = 1.;
	this.XTH_theta                      = 0.;
	this.vertex_equation                = NaN;
	this.element_equation               = NaN;
	this.borderSSA                      = NaN;
	this.borderHO                       = NaN;
	this.borderFS                       = NaN;
	this.setdefaultparameters();
	//}}}
}
