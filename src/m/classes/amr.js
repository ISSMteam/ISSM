//AMR class definition
//
//   Usage:
//      amr=new amr();

function amr (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//hmin and hmax
		this.hmin								= 100.;
		this.hmax								= 100.e3;

		//fields
		this.fieldname							= "Vel";
		this.err								= 3.;

		//keep metric?
		this.keepmetric							= 1;

		//control of element lengths
		this.gradation							= 1.5;

		//other criterias
		this.groundingline_resolution			= 500.;
		this.groundingline_distance				= 0;
		this.icefront_resolution				= 500;
		this.icefront_distance					= 0;
		this.thicknesserror_resolution			= 500;
		this.thicknesserror_threshold			= 0;
		this.thicknesserror_groupthreshold 		= 0;
		this.thicknesserror_maximum				= 0;
		this.deviatoricerror_resolution			= 500;	
		this.deviatoricerror_threshold			= 0;	
		this.deviatoricerror_groupthreshold		= 0;	
		this.deviatoricerror_maximum			= 0;	

		//is restart? This calls femmodel->ReMesh() before first time step. 
		this.restart							= 0;
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   amr parameters:'));
		fielddisplay(this,'hmin','minimum element length');
		fielddisplay(this,'hmax','maximum element length');
		fielddisplay(this,'fieldname','name of input that will be used to compute the metric (should be an input of FemModel)');
		fielddisplay(this,'keepmetric','indicates whether the metric should be kept every remeshing time');
		fielddisplay(this,'gradation','maximum ratio between two adjacent edges');
		fielddisplay(this,'groundingline_resolution','element length near the grounding line');
		fielddisplay(this,'groundingline_distance','distance around the grounding line which elements will be refined');
		fielddisplay(this,'icefront_resolution','element length near the ice front');
		fielddisplay(this,'icefront_distance','distance around the ice front which elements will be refined');
		fielddisplay(this,'thicknesserror_resolution','element length when thickness error estimator is used');
		fielddisplay(this,'thicknesserror_threshold','maximum threshold thickness error permitted');
		fielddisplay(this,'thicknesserror_groupthreshold','maximum group threshold thickness error permitted');
		fielddisplay(this,'thicknesserror_maximum','maximum thickness error permitted');
		fielddisplay(this,'deviatoricerror_resolution','element length when deviatoric stress error estimator is used');
		fielddisplay(this,'deviatoricerror_threshold','maximum threshold deviatoricstress error permitted');
		fielddisplay(this,'deviatoricerror_groupthreshold','maximum group threshold deviatoricstress error permitted');
		fielddisplay(this,'deviatoricerror_maximum','maximum deviatoricstress error permitted');
		fielddisplay(this,'deviatoricerror_maximum','maximum deviatoricstress error permitted');
		fielddisplay(this,'restart','indicates if ReMesh() will call before first time step');
	}// }}}
	this.classname= function(){// {{{
		return "amr";

	}// }}}
	this.checkconsistency = function(md,solution,analyses) { //{{{
		checkfield(md,'fieldname','amr.hmax','numel',[1],'>',0,'NaN',1);
		checkfield(md,'fieldname','amr.hmin','numel',[1],'>',0,'<',this.hmax,'NaN',1);
		checkfield(md,'fieldname','amr.keepmetric','numel',[1],'>=',0,'<=',1,'NaN',1);
		checkfield(md,'fieldname','amr.gradation','numel',[1],'>=',1.1,'<=',5,'NaN',1);
		checkfield(md,'fieldname','amr.groundingline_resolution','numel',[1],'>',0,'<',this.hmax,'NaN',1);
		checkfield(md,'fieldname','amr.groundingline_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','amr.icefront_resolution','numel',[1],'>',0,'<',this.hmax,'NaN',1);
		checkfield(md,'fieldname','amr.icefront_distance','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','amr.thicknesserror_resolution','numel',[1],'>',0,'<',this.hmax,'NaN',1);
		checkfield(md,'fieldname','amr.thicknesserror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);
		checkfield(md,'fieldname','amr.thicknesserror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);
		checkfield(md,'fieldname','amr.thicknesserror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','amr.deviatoricerror_resolution','numel',[1],'>',0,'<',this.hmax,'NaN',1);
		checkfield(md,'fieldname','amr.deviatoricerror_threshold','numel',[1],'>=',0,'<=',1,'NaN',1);
		checkfield(md,'fieldname','amr.deviatoricerror_groupthreshold','numel',[1],'>=',0,'<=',1,'NaN',1);
		checkfield(md,'fieldname','amr.deviatoricerror_maximum','numel',[1],'>=',0,'NaN',1,'Inf',1);
		checkfield(md,'fieldname','amr.restart','numel',[1],'>=',0,'<=',1,'NaN',1);
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{
		WriteData(fid,prefix,'name','md.amr.type','data',1,'format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','hmin','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','hmax','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','fieldname','format','String');
		WriteData(fid,prefix,'object',this,'fieldname','err','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','keepmetric','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','gradation','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','groundingline_resolution','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','groundingline_distance','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','icefront_resolution','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','icefront_distance','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','thicknesserror_resolution','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','thicknesserror_threshold','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','thicknesserror_groupthreshold','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','thicknesserror_maximum','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','deviatoricerror_resolution','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','deviatoricerror_threshold','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','deviatoricerror_groupthreshold','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','deviatoricerror_maximum','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','restart','format','Integer');
	}//}}}
	this.fix=function() { //{{{
	}//}}}
	//properties 
	// {{{
	this.hmin								= 0.;
	this.hmax								= 0.;
	this.fieldname							= "";
	this.err								= 0.;
	this.keepmetric							= 0;
	this.gradation							= 0.;
	this.groundingline_resolution			= 0.;
	this.groundingline_distance				= 0.;
	this.icefront_resolution				= 0.;
	this.icefront_distance					= 0.;
	this.thicknesserror_resolution			= 0.;
	this.thicknesserror_threshold			= 0.;
	this.thicknesserror_groupthreshold		= 0.;
	this.thicknesserror_maximum				= 0.;
	this.deviatoricerror_resolution			= 0.;
	this.deviatoricerror_threshold			= 0.;
	this.deviatoricerror_groupthreshold		= 0.;
	this.deviatoricerror_maximum			= 0.;
	this.restart							= 0.;

	this.setdefaultparameters();
	//}}}
}
