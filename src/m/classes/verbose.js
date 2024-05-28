//VERBOSE class definition
//
//   Available verbosity levels:
//      mprocessor  : model processing 
//      module      : modules
//      solution    : solution sequence
//      solver      : solver info (extensive)
//      convergence : convergence criteria
//      control     : control method
//      qmu         : sensitivity analysis
//      autodiff    : AD analysis
//      smb         : smb analysis
//
//   Usage:
//      verbose=verbose();
//      verbose=verbose(3);
//      verbose=verbose('all');
//      verbose=verbose('001100');
//      verbose=verbose('module',true,'solver',false);

//WARNING: some parts of this file are Synchronized with src/c/shared/Numerics/Verbosity.h
//         Do not modify these sections. See src/c/shared/Numerics/README for more info

function verbose (){
	//methods
	this.setdefaultparameters = function(){// {{{
		//switch(nargin),
			//case 0,
		this.verbose.solution=true;
		this.verbose.qmu=true;
		this.verbose.control=true;
	}// }}}
	this.disp= function(){// {{{
		//BEGINDISP
		console.log(sprintf('verbose class echo:'));
		console.log(sprintf('   %s : %i','mprocessor',this.mprocessor));
		console.log(sprintf('   %s : %i','module',this.module));
		console.log(sprintf('   %s : %i','solution',this.solution));
		console.log(sprintf('   %s : %i','solver',this.solver));
		console.log(sprintf('   %s : %i','convergence',this.convergence));
		console.log(sprintf('   %s : %i','control',this.control));
		console.log(sprintf('   %s : %i','qmu',this.qmu));
		console.log(sprintf('   %s : %i','autodiff',this.autodiff));
		console.log(sprintf('   %s : %i','smb',this.smb));
		//ENDDISP
	}// }}}
		this.checkconsistency = function(md,solution,analyses){ // {{{

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'data',this.VerboseToBinary(),'name','md.verbose','format','Integer');
		}//}}}
		this.VerboseToBinary = function () { //{{{

			//BEGINVERB2BIN
			var binary=0;
			if (this.mprocessor) binary=binary|1; 
			if (this.module) binary=binary|2; 
			if (this.solution) binary=binary|4; 
			if (this.solver) binary=binary|8; 
			if (this.convergence) binary=binary|16; 
			if (this.control) binary=binary|32; 
			if (this.qmu) binary=binary|64; 
			if (this.autodiff) binary=binary|128; 
			if (this.smb) binary=binary|256; 
			//ENDVERB2BIN
			return binary;

		} //}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.mprocessor=false;
	this.module=false;
	this.solution=false;
	this.solver=false;
	this.convergence=false;
	this.control=false;
	this.qmu=false;
	this.autodiff=false;
	this.smb=false;

	this.setdefaultparameters();
	// }}}
}
