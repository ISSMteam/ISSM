//LOCAL cluster class definition
//
//   Usage:
//      local=new local();

function local (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   local cluster class echo: []'));
	}// }}}
	this.classname= function(){// {{{
		return "local";
	}// }}}
		this.checkconsistency = function (md,solution,analyses) { //{{{
		} //}}}
	//properties 
	// {{{
	this.setdefaultparameters();
	//}}}
}
