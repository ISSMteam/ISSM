//TOOLKITS class definition
//
//   Usage:
//      toolkits=new toolkits();

function toolkits (){
	//methods
	this.setdefaultparameters = function(){// {{{

		//default toolkits: 
		if (IssmConfig('_HAVE_PETSC_')){
			//MUMPS is the default toolkits
			if (IssmConfig('_HAVE_MUMPS_')){
				this.DefaultAnalysis           = mumpsoptions();
			}
			else{
				this.DefaultAnalysis           = iluasmoptions(); 
			}
		}
		else{
			if (IssmConfig('_HAVE_MUMPS_')){
				this.DefaultAnalysis           = issmmumpssolver(); 
			}
			else if (IssmConfig('_HAVE_GSL_')){
				this.DefaultAnalysis           = issmgslsolver(); 
			}
			else{
				console.warn('toolkits setdefaultparameters message: need at least Mumps or Gsl to define an issm solver type, no default solver assigned');
			}
		}

		this.RecoveryAnalysis = this.DefaultAnalysis;
	}// }}}
	this.disp = function(){// {{{
		console.log(sprintf('List of toolkits options per analysis:\n'));
		for(var prop in this){
			if(typeof this[prop] == 'object'){
				console.log(prop+ ':',this[prop]);
			}
		}
	}// }}}
	this.checkconsistency = function (md,solution,analyses) { // {{{
		for(var prop in this){
			if(typeof this[prop] == 'object'){
				if (this[prop] == ''){
					md.checkmessage(sprintf("md.toolkits.%s is empty",prop));
				}
			}
		}
	} // }}}
		 this.ToolkitsFile = function(filename) { //{{{
		 //TOOLKITSFILE - build toolkits file (in string format)
		 //
		 //   Build a Petsc compatible options string, from the toolkits model field  + return options string. 
		 //   This file string will also be used when the toolkit used is 'issm' instead of 'petsc'
		 //
		 //   Usage:     var toolkitsstring = toolkits.ToolkitsFile();

			 var string = '';

			 //write header
			 string += sprintf('%s%s%s\n','\%Toolkits options file: ',filename,' written from Javascript toolkits array');

			 //start writing options
			 for (var analysis in this){
				 var options;
				 
				 if(typeof this[analysis] == 'object') options=this[analysis]; else continue;

				 //first write analysis:
				 string += sprintf('\n+%s\n',analysis); //append a + to recognize it's an analysis string

				 //now, write options
			
				 for(var optionname in options){
					 var optionvalue=options[optionname];

					 if (optionvalue.length==0){
						 //this option has only one argument
						 string+=sprintf('-%s\n',optionname);
					 }
					 else{
						 //option with value. value can be string or scalar
						 if (typeof optionvalue == 'number'){
							 string+=sprintf('-%s %g\n',optionname,optionvalue);
						 }
						 else if (typeof optionvalue == 'string'){
							 string+=sprintf('-%s %s\n',optionname,optionvalue);
						 }
						 else throw Error(sprintf("ToolkitsFile error: option '%s' is not well formatted",optionname));
					 }
				 }
			 }
			 return string;
		 } //}}}
	//properties 
	// {{{
	this.DefaultAnalysis  = [];
	this.RecoveryAnalysis = [];
	//The other properties are dynamic
	this.setdefaultparameters();
	//}}}
}
