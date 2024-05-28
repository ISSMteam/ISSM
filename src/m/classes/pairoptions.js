//PAIROPTIONS class definition
//
//   Usage:
//      pairoptions=pairoptions();
//      pairoptions=pairoptions('module',true,'solver',false);

function pairoptions(args) { 


	//methods
	this.constructor = function (args) {// {{{

		//initialize list
		if (args.length==0){
			this.list=[];
		}
		else{

			//check length of input
			if (args.length % 2){
				throw Error('pairoptions error message: invalid parameter/value pair arguments') 
			}
			numoptions = args.length/2;

			//Allocate memory
			this.list= Create2DArray(numoptions,3);

			//go through args and build list of obj
			for (var i=0;i<numoptions;i++){
				if (typeof args[2*i] === 'string'){
					this.list[i][0]=args[2*i];
					this.list[i][1]=args[2*i+1];
					this.list[i][2]=false;  //used?
				}
				else{
					//option is not a string, ignore it
					console.log(sprintf('%s%i%s\n','WARNING: option number ',i,' is not a string, it will be ignored'));
					this.list[i][0]=[];
					this.list[i][1]=[];
					this.list[i][2]=[];
					continue
				}
			}
		}
	}// }}}
	this.addfield = function (field, value){ // {{{
		if (typeof field === 'string'){
			this.list.push([field,value,false]);
		}
	}// }}}
	this.numoptions = function (){ // {{{
		return this.list.length;
	}// }}}
	this.addfielddefault = function (field,value){ // {{{
		//ADDFIELDDEFAULT - add a field to an options list if it does not exist
		if (typeof field === 'string'){
			if (!this.exist(field)){
				this.list.push([field,value,true]); //true is a default so user will not be notified if not used
			}
		}
	} // }}}
	this.AssignObjectFields = function(object){ // {{{
		//ASSIGNOBJECTFIELDS - assign object fields from options
		for (var i=0;i<list.length;i++){
			fieldname=list[i][0];
			fieldvalue=list[i][1];
			if (fieldname in object){
				obj2[fieldname]=fieldvalue;
			}
			else{
				console.log(sprintf("%s'%s'%s%s\n",'WARNING: ',fieldname, 'is not a property of ',typeof object));
			}
		}
	} // }}}
	this.changefieldvalue = function(field,newvalue){ // {{{
		//CHANGEOPTIONVALUE - change the value of an option in an option list

		var found=0;
		for (var i=0;i<this.list.length;i++){
			if (this.list[i][0] === field){
				found=1;
			}
		}

		if (found==0){
			this.list.push([field,newvalue,true]); // do not notify user if unused
		}
		else{
			for (var i=0;i<this.list.length;i++){
				if (this.list[i][0] === field){
					this.list[i][1] = newvalue;
				}
			}
		}
	} // }}}
	this.deleteduplicates = function(warn){ // {{{
		//DELETEDUPLICATES - delete duplicates in an option list

		//track the first occurrence of each option
		var indices=NewArrayFill(this.list.length,0);
		for (var i=0;i<this.list.length;i++){
			if(indices[i]==0){
				for(var j=i+1;j<this.list.length;j++){
					if (this.list[i][0] === this.list[j][0])indices[j]=1;
				}
			}
		}
		sumindices=ArraySum(indices);

		//remove duplicates from the options list
		newlist=Create2DArray(sumindices,3);
		var count=0;
		for (var i=0;i<this.list.length;i++){
			if (indices[i]==1) if (warn) console.log(sprintf("%s%s%s\n",'WARNING: option ', this.list[i,0],' appeared more than once. Only its first occurrence will be kept'));
			else{
				newlist[count]=this.list[i];
				count++;
			}
		}
	} // }}}
	this.displayunused = function (){ // {{{
		//DISPLAYUNUSED - display unused options

		for (var i=0;i<this.list.length;i++){
			if (!(this.list[i][2])){
				console.log(sprintf("%s%s%s\n",'WARNING: option ',this.list[i][0],' was not used'));
			}
		}
	}// }}}
	this.disp = function (){ //{{{
		if (this.list.length){
			console.log(sprintf('   pairoptions: (%i)\n',this.list.length));
			for (var i=0;i<this.list.length;i++){
				if (typeof this.list[i][1] === 'string'){
					console.log(sprintf("     field: '%s' value(string): ''%s''",this.list[i][0],this.list[i][1]));
				}
				else if( typeof this.list[i][1] === 'number'){
					console.log(sprintf("     field: '%s' value(number): %g",this.list[i][0],this.list[i][1]));
				}
				else if( IsArray(this.list[i][1])){
					console.log(sprintf("     field: '%s' value(array): [%i]",this.list[i][0],this.list[i][1].length));
				}
			}
		}
		else{
			console.log(sprintf('   list: empty'));
		}
	}// }}}
	this.exist = function (field) { //{{{

		//EXIST - check if the option exists
		//some argument checking: 
		if (!(typeof field === 'string')){
			throw Error('exist error message: field should be a string');
		}

		//Recover option
		var bool=0;
		for (var i=0;i<this.list.length;i++){
			if (this.list[i][0] === field){
				bool=1;
				this.list[i][2]=1; //It is a default so user will not be notified if not used
				break;
			}
		}
		return bool;
	} // }}}
	this.fieldoccurrences = function(field){ // {{{

		//FIELDOCCURRENCES - get number of occurrence of a field
		var num=0;

		//check input 
		if (!(typeof field === 'string')){
			throw Error('exist error message: field should be a string');
		}

		//count number of occurrences:
		for (var i=0;i<this.list.length;i++) if (this.list[i][0] === field)num++;

		return num;

	} // }}}
	this.getfieldvalue = function(field){ // {{{
		//GETOPTION - get the value of an option
		//
		//   Usage:
		//      value=pairoptions.getfieldvalue(field,varargin)
		//
		//   Find an option value from a field. A default option
		//   can be given in input if the field does not exist
		//
		//   Examples:
		//      value=pairoptions.getfieldvalue('caxis');
		//      value=pairoptions.getfieldvalue('caxis',[0 2]);

		//some argument checking: 
		if(!(arguments.length==1 | arguments.length==2)){
			error('pairoptions usage error: getfieldvalue bad usage');
		}

		if (!(typeof field === 'string')){
			throw Error('pairoptions error message: field should be a string');
		}

		//Recover option
		for(var i=0;i<this.list.length;i++){
			if (this.list[i][0] === field){
				this.list[i][2]=1; //option used
				return value=this.list[i][1];
			}
		}

		//The option has not been found, output default if provided
		if (arguments.length==2){
			return arguments[1];
		}
		else{
			throw Error(sprintf("%s%s%s\n",'error message: field ',field,' has not been provided by user (and no default value has been specified)'));
		}
	} // }}}
	this.removefield = function(field,warn){// {{{

		//REMOVEFIELD - delete a field in an option list
		//
		//   Usage:
		//      options.removefield(field,warn)
		//
		//   if warn==1 display an info message to warn user that
		//   some of his options have been removed.

		//check if field exists
		if (this.exist(field)){

			var indices;
			var count;

			//find where the field is located
			indices=NewArrayFill(this.list.length,1);
			for (var i=0;i<this.list.length;i++)if(this.list[i][1] === field)indices[i]=0;
			sumindices=ArraySum(indices);

			//remove duplicates from the options list
			newlist=Create2DArray(sumindices,3);

			count=0;
			for (var i=0;i<this.list.length;i++){
				if(!(this.list[i][1] === field)){
					newlist[count]=this.list[i];
					count++;
				}
			}
			this.list=newlist;

			//warn user if requested
			if (warn){
				console.log(sprintf("%s%s%s\n",'removefield info: option ',field,' has been removed from the list of options.'));
			}
		}
	} // }}}
	this.marshall = function(fid,firstindex){// {{{

		throw Error('pairoptions marshall error: not implemented yet!');
	} // }}}

	//properties 
	this.list         = [];
	this.constructor(args);
}
