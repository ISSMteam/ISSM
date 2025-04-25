//AUTODIFF class definition
//
//   Usage:
//      autodiff=new autodiff();

function autodiff (){
	//methods
	this.setdefaultparameters = function(){// {{{

		this.obufsize     = 524288;
		this.lbufsize     = 524288;
		this.cbufsize     = 524288;
		this.tbufsize     = 524288;
		this.gcTriggerRatio=2.0;
		this.gcTriggerMaxSize=65536;
		this.tapeAlloc    = 15000000;
		this.outputTapeMemory = 0;
		this.outputTime = 0;
		this.enablePreaccumulation = 0;

	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   automatic differentiation parameters:'));
		fielddisplay(this,'isautodiff','indicates if the automatic differentiation is activated');
		fielddisplay(this,'dependents','list of dependent variables');
		fielddisplay(this,'independents','list of independent variables');
		fielddisplay(this,'driver',"ADOLC driver ('fos_forward' or 'fov_forward')");
		fielddisplay(this,'obufsize','Number of operations per buffer (==OBUFSIZE in usrparms.h)');
		fielddisplay(this,'lbufsize','Number of locations per buffer (==LBUFSIZE in usrparms.h)');
		fielddisplay(this,'cbufsize','Number of values per buffer (==CBUFSIZE in usrparms.h)');
		fielddisplay(this,'tbufsize','Number of taylors per buffer (<=TBUFSIZE in usrparms.h)');
		fielddisplay(this,'gcTriggerRatio','free location block sorting/consolidation triggered if the ratio between allocated and used locations exceeds gcTriggerRatio');
		fielddisplay(this,'gcTriggerMaxSize','free location block sorting/consolidation triggered if the allocated locations exceed gcTriggerMaxSize');
		fielddisplay(this,'tapeAlloc','Iteration count of a priori memory allocation of the AD tape');
		fielddisplay(this,'outputTapeMemory','Write AD tape memory statistics to file ad_mem.dat');
		fielddisplay(this,'outputTime','Write AD recording and evaluation times to file ad_time.dat');
		fielddisplay(this,'enablePreaccumulation','Enable CoDiPack preaccumulation in augmented places');

	}// }}}
	this.classname= function(){// {{{
		return "autodiff";
	}// }}}
	this.checkconsistency = function(md,solution,analyses){ //{{{

		//Early return 
		if (!this.isautodiff) return; 

		//Driver value:
		checkfield(md,'fieldname','autodiff.driver','values',['fos_forward','fov_forward','fov_forward_all','fos_reverse','fov_reverse','fov_reverse_all']);

		//buffer values: 
		checkfield(md,'fieldname','autodiff.obufsize','>=',16);
		checkfield(md,'fieldname','autodiff.lbufsize','>=',16);
		checkfield(md,'fieldname','autodiff.cbufsize','>=',16);
		checkfield(md,'fieldname','autodiff.tbufsize','>=',16);
		checkfield(md,'fieldname','autodiff.gcTriggerRatio','>=',0);
		checkfield(md,'fieldname','autodiff.gcTriggerMaxSize','>=',65536);
		checkfield(md,'fieldname','autodiff.tapeAlloc','>=',0);

		//Memory and time output
		checkfield(md,'fieldname','autodiff.outputTapeMemory','numel', [1], 'values', [0, 1]);
		checkfield(md,'fieldname','autodiff.outputTime','numel', [1], 'values', [0, 1]);

		//Memory reduction options
		checkfield(md,'fieldname','autodiff.enablePreaccumulation','>=',0);

		//go through our dependents and independents and check consistency: 
		for (var i=0;i<this.dependents.length;i++){
			dep=this.dependents[i];
			dep.checkconsistency(md,solution,analyses);
		}
		for (var i=0;i<this.independents.length;i++){
			indep=this.independents[i];
			indep.checkconsistency(md,i,solution,analyses,this.driver);
		}
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{

		WriteData(fid,prefix,'object',this,'fieldname','isautodiff','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','driver','format','String');

		//early return
		if (!this.isautodiff){
			WriteData(fid,prefix,'data',false,'name','md.autodiff.mass_flux_segments_present','format','Boolean');
			WriteData(fid,prefix,'data',false,'name','md.autodiff.keep','format','Boolean');
			return;
		}

		//buffer sizes {{{
		WriteData(fid,prefix,'object',this,'fieldname','obufsize','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','lbufsize','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','cbufsize','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','tbufsize','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','gcTriggerRatio','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','gcTriggerMaxSize','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','tapeAlloc','format','Integer');
		//}}}
		//output of memory and time {{{
		WriteData(fid,prefix,'object',this,'fieldname','outputTapeMemory','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','outputTime','format','Boolean');
		//}}}
		//memory reduction options {{{
		WriteData(fid,prefix,'object',this,'fieldname','enablePreaccumulation','format','Boolean');
		//}}}
		//process dependent variables {{{
		num_dependent_objects=this.dependents.length;
		WriteData(fid,prefix,'data',num_dependent_objects,'name','md.autodiff.num_dependent_objects','format','Integer');

		if(num_dependent_objects){
			var names=[];
			types=NewArrayFill(num_dependent_objects,0);
			indices=NewArrayFill(num_dependent_objects,0);

			for (var i=0;i<num_dependent_objects;i++){
				dep=this.dependents[i];

				names.push(dep.name);
				indices[i]=dep.index;
			}
			WriteData(fid,prefix,'data',names,'name','md.autodiff.dependent_object_names','format','StringArray');
			WriteData(fid,prefix,'data',indices,'name','md.autodiff.dependent_object_indices','format','IntMat','mattype',3);
		}
		//}}}
		//process independent variables {{{
		num_independent_objects=this.independents.length;
		WriteData(fid,prefix,'data',num_independent_objects,'name','md.autodiff.num_independent_objects','format','Integer');

		if(num_independent_objects){
			names=NewArrayFill(num_independent_objects,0);
			types=NewArrayFill(num_independent_objects,0);

			for (var i=0;i<num_independent_objects;i++){
				indep=this.independents[i];

				names[i]=indep.name;
				types[i]=indep.typetoscalar();
			}
			WriteData(fid,prefix,'data',names,'name','md.autodiff.independent_object_names','format','StringArray');
			WriteData(fid,prefix,'data',types,'name','md.autodiff.independent_object_types','format','IntMat','mattype',3);
		}
		//}}}
		//if driver is fos_forward, build index:  {{{
		if (this.driver == 'fos_forward'){
			var index=0;

			for (var i=0;i<num_independent_objects;i++){
				indep=this.independents[i];
				if (!(isNaN(indep.fos_forward_index))){
					index=index+indep.fos_forward_index;
					break;
				}
				else{
					if (indep.type=='scalar') index=index+1;
					else index=index+indep.nods;
				}
			}
			index=index-1; //get c-index numbering going
			WriteData(fid,prefix,'data',index,'name','md.autodiff.fos_forward_index','format','Integer');
		}
		//}}}
		//if driver is fos_reverse, build index:  {{{
		if (this.driver  == 'fos_reverse'){
			var index=0;

			for (var i=0;i<num_dependent_objects;i++){
				dep=this.dependents[i];
				if (!(isNaN(dep.fos_reverse_index))){
					index=index+dep.fos_reverse_index;
					break;
				}
				else{
					if (dep.type =='scalar') index=index+1;
					else index=index+dep.nods;
				}
			}
			index=index-1; //get c-index numbering going
			WriteData(fid,prefix,'data',index,'name','md.autodiff.fos_reverse_index','format','Integer');
		}
		//}}}
		//if driver is fov_forward, build indices:  {{{
		if (this.driver == 'fov_forward'){
			var indices=0;

			for (var i=0;i<num_independent_objects;i++){
				indep=this.independents[i];
				if (!indep.fos_forward_index.length){
					indices=indices+indep.fov_forward_indices;
					break;
				}
				else{
					if (indep.type =='scalar') indices=indices+1;
					else indices=indices+indep.nods;
				}
			}
			indices=indices-1; //get c-indices numbering going
			WriteData(fid,prefix,'data',indices,'name','md.autodiff.fov_forward_indices','format','IntMat','mattype',3);
		}
		//}}}
		//deal with mass fluxes:  {{{
		mass_flux_segments=[];
		for (var i=0;i<num_dependent_objects;i++){
			dep=this.dependents[i];
			if (dep.name =='MassFlux'){
				mass_flux_segments.push(dep.segments);
			}
		}
		if (mass_flux_segments.length){
			WriteData(fid,prefix,'data',mass_flux_segments,'name','md.autodiff.mass_flux_segments','format','MatArray');
			flag=true;
		}
		else flag=false;
		WriteData(fid,prefix,'data',flag,'name','md.autodiff.mass_flux_segments_present','format','Boolean');
		//}}}
		//deal with trace keep on: {{{
		keep=false;

		//From ADOLC userdoc: 
		// The optional integer argument keep of trace on determines whether the numerical values of all active variables are 
		// recorded in a buffered temporary array or file called the taylor stack. This option takes effect if keep = 1 and 
		// prepares the scene for an immediately following gradient evaluation by a call to a routine implementing the reverse 
		// mode as described in the Section 4 and Section 5. 
		//

		if (this.driver.length<=3) keep=false; //there is no "_reverse" string within the driver string: 
		else{
			if (this.driver.splice(4) == '_reverse') keep=true;
			else keep=false;
		}
		WriteData(fid,prefix,'data',keep,'name','md.autodiff.keep','format','Boolean');
		//}}}
	}//}}}
	this.fix=function() { //{{{
		this.obufsize=NullFix(this.obufsize,NaN);
		this.lbufsize=NullFix(this.lbufsize,NaN);
		this.cbufsize=NullFix(this.cbufsize,NaN);
		this.tbufsize=NullFix(this.tbufsize,NaN);
		this.gcTriggerRatio=NullFix(this.gcTriggerRatio,NaN);
		this.gcTriggerMaxSize=NullFix(this.gcTriggerMaxSize,NaN);
		this.tapeAlloc=NullFix(this.tapeAlloc,NaN);
	}//}}}
	//properties 
	// {{{
	this.isautodiff   = false;
	this.dependents   = [];
	this.independents = [];
	this.driver       = 'fos_forward';
	this.obufsize     = NaN;
	this.lbufsize     = NaN;
	this.cbufsize     = NaN;
	this.tbufsize     = NaN;
	this.gcTriggerRatio = NaN;
	this.gcTriggerMaxSize = NaN;
	this.tapeAlloc = NaN;

	this.setdefaultparameters();
	//}}}
}
