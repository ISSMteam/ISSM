//OUTPUTDEFINITION class definition
//
//   Usage:
//      outputdefinition=new outputdefinition();

function outputdefinition (){
	//methods
	this.setdefaultparameters = function(){// {{{
		this.definitions=[];
	}// }}}
	this.disp= function(){// {{{
		console.log(sprintf('   outputdefinition:'));
		fielddisplay(this,'definitions','list of potential outputs that can be requested, but which need additional data to be defined');


	}// }}}
	this.classname= function(){// {{{
		return "outputdefinition";
	}// }}}
		this.checkconsistency = function(md,solution,analyses) { //{{{

			checkfield(md,'fieldname','outputdefinition.definitions','cell',1);

			for (var i=0;i<this.definitions.length;i++){
				this.definitions[i].checkconsistency(md,solution,analyses);
			}

		} // }}}
		this.marshall=function(md,prefix,fid) { //{{{

			var data=NewArrayFill(this.definitions.length,'');	
			for(var i=0;i<this.definitions.length;i++){
				this.definitions[i].marshall(md,fid,prefix);
				classdefinition=this.definitions[i].classname();
				classdefinition=classdefinition.charAt(0).toUpperCase() + classdefinition.slice(1); // we match our string definitions
				data[i]=classdefinition;
			}
			data=ArrayUnique(data);
			if(data.length==0){ data=''; }

			WriteData(fid,prefix,'data',data,'name','md.outputdefinition.list','format','StringArray');
		}//}}}
		this.fix=function() { //{{{
		}//}}}
	//properties 
	// {{{
	this.definitions                 = [];
	this.setdefaultparameters();
	//}}}
}
