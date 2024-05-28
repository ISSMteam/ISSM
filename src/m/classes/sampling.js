//SAMPLING class definition
//
//   Usage:
//      sampling=new sampling();

function sampling (){
	//methods
	this.defaultoutputs = function(){// {{{
       	return [];
	}// }}}
	this.setdefaultparameters = function(){// {{{

		//Scaling coefficient
		this.tau=1;  

		//Apply Robin boundary conditions
		this.robin=0;   
		
		//Temporal correlation factor
		this.phi=0;  
		
		//Exponent in fraction SPDE (default=2, biLaplacian covariance
		//operator)
		this.alpha=2; // Default 
		
		//Seed for pseudorandom number generator (default -1 for random seed)
		this.seed=-1;
		
		//default output
		this.requested_outputs=['default'];
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   Sampling parameters:'));

		console.log(sprintf('\n      %s','Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):'));
		fielddisplay(this,'kappa','coefficient of the identity operator');
		fielddisplay(this,'tau','scaling coefficient of the solution (default 1.0)');
		fielddisplay(this,'alpha','exponent in PDE operator, (default 2.0, BiLaplacian covariance operator)');
	  
		console.log(sprintf('\n      %s','Parameters of Robin boundary conditions nabla () \cdot normvec + beta ():'));
		fielddisplay(this,'robin','Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default 0)');
		fielddisplay(this,'beta','Coefficient in Robin boundary conditions (to be defined for robin = 1)');          
		
		console.log(sprintf('\n      %s','Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):'));
		fielddisplay(this,'phi','Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)');
		
		console.log(sprintf('\n      %s','Other parameters of stochastic sampler:'));
		fielddisplay(this,'seed','Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default -1)');
		fielddisplay(this,'requested_outputs','additional outputs requested (not implemented yet)');
	}// }}}
	this.classname= function(){// {{{
		return "sampling";
	}// }}}
	this.checkconsistency = function(md,solution,analyses) { //{{{
		if (!ismember('SamplingAnalysis',analyses)) return;
		
		checkfield(md,'fieldname','sampling.kappa','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0);
		checkfield(md,'fieldname','sampling.tau','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','sampling.robin','numel',1,'values',[0,1]);
		if (md.sampling.robin) {
			checkfield(md,'fieldname','sampling.beta','NaN',1,'Inf',1,'size',[md.mesh.numberofvertices,1],'>',0);
		}
		checkfield(md,'fieldname','sampling.phi','NaN',1,'Inf',1,'numel',1,'>=',0);
		checkfield(md,'fieldname','sampling.alpha','NaN',1,'Inf',1,'numel',1,'>',0);
		checkfield(md,'fieldname','sampling.seed','NaN',1,'Inf',1,'numel',1);
		checkfield(md,'fieldname','sampling.requested_outputs','stringrow',1);
	} // }}}
	this.marshall=function(md,prefix,fid) { //{{{
		WriteData(fid,prefix,'object',this,'fieldname','kappa','format','DoubleMat','mattype',1);
		WriteData(fid,prefix,'object',this,'fieldname','tau','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','beta','format','DoubleMat','mattype',1);
		WriteData(fid,prefix,'object',this,'fieldname','phi','format','Double');
		WriteData(fid,prefix,'object',this,'fieldname','alpha','format','Integer');
		WriteData(fid,prefix,'object',this,'fieldname','robin','format','Boolean');
		WriteData(fid,prefix,'object',this,'fieldname','seed','format','Integer');
		
		//process requested outputs
		outputs = this.requested_outputs;
		pos = find(ismember(outputs,'default'));
		if (!isempty(pos)) {
			ArrayIndex(outputs,pos,[]);                         //remove 'default' from outputs
			outputs = ArrayConcat(outputs, this.defaultoutputs()); //add defaults
		}
		WriteData(fid,prefix,'data',outputs,'name','md.sampling.requested_outputs','format','StringArray');
	}//}}}
	//properties 
	// {{{
    this.kappa          	= NaN;
	this.tau            	= 0;
	this.beta           	= NaN;
    this.phi            	= 0;
    this.alpha          	= 0;
    this.robin          	= 0;
    this.seed           	= 0;
	this.requested_outputs  = [];

	this.setdefaultparameters();
	//}}}
}
