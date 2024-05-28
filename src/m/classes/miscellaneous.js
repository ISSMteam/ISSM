//MISCELLANEOUS class definition
//
//   Usage:
//      miscellaneous=new miscellaneous();

function miscellaneous (){
	//methods
	this.setdefaultparameters = function(){// {{{
	}// }}}
	this.disp= function(){// {{{

		console.log(sprintf('   miscellaneous parameters:'));

		fielddisplay(this,'notes','notes in a cell of strings');
		fielddisplay(this,'name','model name');
		fielddisplay(this,'dummy','empty field to store some data');

	}// }}}
	this.classname= function(){// {{{
		return "miscellaneous";
	}// }}}
		this.checkconsistency= function(md,solution,analyses) {// {{{

			checkfield(md,'fieldname','miscellaneous.name','empty',1);
		}// }}}
		this.marshall=function(md,prefix,fid) { //{{{
			WriteData(fid,prefix,'object',this,'fieldname','name','format','String');
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.notes = '';
	this.name  = '';
	this.dummy = [];

	this.setdefaultparameters();
	//}}}
}
