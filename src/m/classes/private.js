//PRIV class definition
//
//   Usage:
//      priv =new priv();

function priv (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   private parameters: do not change'));

		fielddisplay(this,'isconsistent','is model this consistent');
		fielddisplay(this,'runtimename','name of the run launched');
		fielddisplay(this,'bamg','structure with mesh properties constructed if bamg is used to mesh the domain');
		fielddisplay(this,'solution','type of solution launched');
	}// }}}
	this.checkconsistency = function(md,solution,analyses){ // {{{

	}// % }}}
	//properties 
	// {{{
	this.isconsistent = true;
	this.runtimename  = '';
	this.bamg         = {};
	this.solution     = '';

	this.setdefaultparameters();
	//}}}
}
