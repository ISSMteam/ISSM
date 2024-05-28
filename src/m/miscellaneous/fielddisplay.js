function fielddisplay(md,name,comment){
//FIELDDISPLAY - display model field
//
//   Usage:
//      fielddisplay(md,name,comment)

	//get field
	field=md[name];

	//disp corresponding line as a function of field type (offset set as 9 spaces)
	parsedisplay('         ',name,field,comment);
}

function parsedisplay(offset,name,field,comment) { //{{{

	//string
	if (typeof(field) == "string"){

		if (field.length > 30){
			displayunit(offset,name,"not displayed",comment);
		}
		else{
			displayunit(offset,name,"'"+field+"'",comment);
		}
	}
	//numeric
	else if (typeof(field) == "number"){
		
		displayunit(offset,name,sprintf("%g",field),comment);

	}
	//logical
	else if (typeof(field) == "boolean") {

		if (field){
			displayunit(offset,name,"true",comment);
		}
		else{
			displayunit(offset,name,"false",comment);
		}

	}
	//object
	else if (typeof(field) == "object"){

		if(field.length == 0) displayunit(offset,name,sprintf("(%i)",field.length),comment);
		else if ((field[0].length==0) | (typeof field[0].length =='undefined')){
			displayunit(offset,name,sprintf("(%i)",field.length),comment);
		}
		else{
			displayunit(offset,name,sprintf("(%i,%i)",field.length,field[0].length),comment);
		}

	}
	else{
		displayunit(offset,name,"not displayed",comment);
	}
} //}}}

function displayunit(offset,name,characterization,comment){ // {{{

	//take care of name
	if (name.length>23){
		name=name.slice(0,21) + "...";
	}

	//take care of characterization
	if ( characterization == "\" \"" || characterization == "NaN" ){
	
		characterization="N/A";
	}
	if (characterization.length>15){
		characterization=characterization.slice(0,13) + "...";
	}

	//print
	if (comment.length==0){
		console.log(sprintf("%s%-23s: %-15s",offset,name,characterization));
	}
	else{
		if (typeof(comment) == "string"){
			//console.log(sprintf("%s%-23s: %-15s -- %s",offset,name,characterization,comment));
			console.log(sprintf("%s%s: %-15s -- %s",offset,name,characterization,comment));
		}
		else{
			throw Error("fielddisplay error message: format for comment not supported yet");
		}
	}
} //}}}
