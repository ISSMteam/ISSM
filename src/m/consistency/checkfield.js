function checkfield(md){
//CHECKFIELD - check field consistency
//
//   Used to check model consistency.
//   Requires: 
//     'field' or 'fieldname' option. If 'fieldname' is provided, it will retrieve it from the model md. (md.(fieldname)) 
//             If 'field' is provided, it will assume the argument following 'field' is a numeric array.
//   Available options:
//      - NaN: 1 if check that there is no NaN
//      - Inf: 1 if check that there is no Inf
//      - size: [lines cols], NaN for non checked dimensions
//      - >:  greater than provided value
//      - >=: greater or equal to provided value
//      - <:  smallerthan provided value
//      - <=: smaller or equal to provided value
//      - < vec:  smallerthan provided values on each vertex
//      - timeseries: 1 if check time series consistency (size and time)
//      - values: array of strings or vector of acceptable values
//      - numel: list of acceptable number of elements
//      - array: 1 if check that is array
//      - empty: 1 if check that non empty
//      - message: overloaded error message
//
//   Usage:
//      checkfield(md,fieldname,options);

	//get options
	var args = Array.prototype.slice.call(arguments);
	var  options = new pairoptions(args.slice(1,args.length));
	
	//get field: 
	if (options.exist('field')){
		field=options.getfieldvalue('field'); 
		fieldname=options.getfieldvalue('fieldname','no fieldname'); 
	}
	else{
		fieldname=options.getfieldvalue('fieldname'); 
		eval("field=md." + fieldname + ";");
	}

	//check empty
	if (options.exist('empty')){
		if (field.length == 0){
			md.checkmessage(options.getfieldvalue('message','field ' + "'" + fieldname + "'" + 'is empty'));
		}
	}

	//Check size
	if (options.exist('size')){
		fieldsize=options.getfieldvalue('size');
		if (isNaN(fieldsize[0])){
			if (field[0].length !=fieldsize[1]){
				md.checkmessage(options.getfieldvalue('message', sprintf("field '%s' should have %i columns",fieldname,fieldsize[1])));
			}
		}
		else if (isNaN(fieldsize[1])){
			if (field.length!= fieldsize[0]){
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have %i lines",fieldname,fieldsize[0])));
			}
		}
		else{
			if (IsArray(field)){
				if ((field.length!=fieldsize[0])){
					md.checkmessage(options.getfieldvalue('message', sprintf("field '%s' should be %ix%i",fieldname,fieldsize[0],fieldsize[1])));
				}
			}
			if (IsArray(field[0])){
				if(field[0].length!=fieldsize[1]){
					md.checkmessage(options.getfieldvalue('message', sprintf("field '%s' should be %ix%i",fieldname,fieldsize[0],fieldsize[1])));
				}
			}
		}
	}

	//Check numel
	if (options.exist('numel')){
		fieldnumel=options.getfieldvalue('numel');
		if (!ArrayIsMember(fieldnumel,[field.length])){
			if (fieldnumel.length==1){
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' size should be %i",fieldname,fieldnumel[0])));
			}
			else if (fieldnumel.length==2){
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' size should be %i or %i",fieldname,fieldnumel[0],fieldnumel[1])));
			}
			else{
				var string=''; 
				for (var i=0;i<fieldnumel.length;i++)string=sprintf("%s or %i",string,fieldnumel[i]);
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' size should be %s",fieldname,string)));
			}
		}
	}

	//check NaN
	if (options.getfieldvalue('NaN',0)){
		field2=MatrixToList(field);
		if (ArrayAnyEqual(field2,NaN)){
			md.checkmessage(options.getfieldvalue('message',sprintf("NaN values found in field %s",field)));
		}
	}

	//check Inf
	if (options.getfieldvalue('Inf',0)){
		field2=MatrixToList(field);
		if (ArrayAnyEqual(field2,Infinity)){
			md.checkmessage(options.getfieldvalue('message',sprintf("Inf values found in field %s",field)));
		}
	}

	//check arry
	if (options.getfieldvalue('array',0)){
		if (!IsArray(field)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should be an array!",fieldname)));
		}
	}

	//check values
	if (options.exist('values')){
		fieldvalues=options.getfieldvalue('values');
		if (typeof fieldvalues[0]== 'string'){
			if (typeof field == 'string'){
				if(ArrayAnyEqual(ArrayIsMember([field],fieldvalues),0)){
					if (fieldvalues.length==1){
						md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' value should be %s",fieldname,fieldvalues[0])));
					}
					else if (fieldvalues.length==2){
						md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' values should be %s or %s",fieldname,fieldvalues[0],fieldvalues[1])));
					}
					else{
						var string=''; 
						for (var i=0;i<fieldvalues.length;i++)string=sprintf("%s or %s",string,fieldvalues[i]);
						md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values in %s",fieldname,string)));
					}
				}
			}
			else{
				var string=''; for (var i=0;i<fieldvalues.length;i++)string=sprintf("%s or %s",string,fieldvalues[i]);
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values in %s",fieldname,string)));
			}
		}
		else{
			if (typeof field == 'number') field2=MatrixToList([field]);
			else field2=MatrixToList(field);
			if (typeof field2[0] == 'number'){
				if(ArrayAnyEqual(ArrayIsMember(field2,fieldvalues),0)){
					var string=''; for (var i=0;i<fieldvalues.length;i++)string=sprintf("%s or %g",string,fieldvalues[i]);
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values in %s",fieldname,string)));
				}
			}
			else{
				var string=''; for (var i=0;i<fieldvalues.length;i++)string=sprintf("%s or %g",string,fieldvalues[i]);
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should be a number in %s",fieldname,string)));
			}
		}
	}
	
	//check greater
	if (options.exist('>=')){
		lowerbound=options.getfieldvalue('>=');
		field2=MatrixToList(field);
		if (options.getfieldvalue('timeseries',0)) field2=MatrixToList(ArrayCopy(field).splice(-1,1));

		if (ArrayAnyBelowStrict(field2,lowerbound)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values above %g",fieldname,lowerbound)));
		}
	}
	if (options.exist('>')){
		lowerbound=options.getfieldvalue('>');
		field2=MatrixToList(field);
		if (options.getfieldvalue('timeseries',0)) field2=MatrixToList(ArrayCopy(field).splice(-1,1));

		if (ArrayAnyBelowOrEqual(field2,lowerbound)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values above %g",fieldname,lowerbound)));
		}
	}
	
	//check smaller
	if (options.exist('<=')){
		upperbound=options.getfieldvalue('<=');
		field2=MatrixToList(field);
		if (options.getfieldvalue('timeseries',0)) field2=MatrixToList(ArrayCopy(field).splice(-1,1));

		if (ArrayAnyAboveOrEqual(field2,upperbound)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values below %g",fieldname,upperbound)));
		}
	}
	
	
	if (options.exist('<')){
		upperbound=options.getfieldvalue('<');
		field2=MatrixToList(field);
		if (options.getfieldvalue('timeseries',0)) field2=MatrixToList(ArrayCopy(field).splice(-1,1));
		if (ArrayAnyAboveStrict(field2,upperbound)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have values below %g",fieldname,upperbound)));
		}
	}

	//Check row of stringrow
	if (options.getfieldvalue('stringrow',0)){
		if (IsArray(field[0])){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have only one row",field)));
		}
		if (!IsArray(field)){
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should be an array of string",fieldname)));
		}
		else{
			for(var i=0;i<field.length;i++){
				if (!(typeof field[i] == 'string')){
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' values should be a cell of strings",fieldname)));
				}
			}
		}
	}

	//check file
	if (options.getfieldvalue('file',0)){
		/*if ~exist(field,'file')
			md.checkmessage(['file provided in ''' fieldname ''': ''' field ''' does not exist']);
		end*/
		throw Error("checkfield error message: file checking on javascript not supported yet!");
	}

	//Check forcings (size and times)
	if (options.getfieldvalue('timeseries',0)){
		if (field.length==md.mesh.numberofvertices | field.length==md.mesh.numberofelements){
			if (IsArray(field[0])){
				md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have only one column as there are md.mesh.numberofvertices lines",fieldname)));
			}
		}
		else if ((field.length==md.mesh.numberofvertices+1) | (field.length==md.mesh.numberofelements+1)){
			var times=field[field.length-1]; var sorted_times=ArraySort(times);
			for(var i=0;i<times.length;i++){
				if(times[i] !=sorted_times[i]){
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' columns should be sorted chronologically",fieldname)));
					break;
				}
			}
			var timesm=ArrayCopy(times).splice(0,-1); var timesp=ArrayCopy(times).shift();
			for(var i=0;i<timesm.length;i++){
				if(timesm[i]==timesp[i]){
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' columns must not contain duplicate timesteps",fieldname)));
					break;
				}
			}
		}
		else{
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have md.mesh.numberofvertices or md.mesh.numberofvertices+1 lines",fieldname)));
		}
	}

	//Check single value forcings (size and times)
	if (options.getfieldvalue('singletimeseries',0)){
		if (field.length==2){
			var times=field[1]; var sorted_times=ArraySort(times);
			for(var i=0;i<times.length;i++){
				if(times[i] !=sorted_times[i]){
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' columns should be sorted chronologically",fieldname)));
					break;
				}
			}
			var timesm=ArrayCopy(times).splice(0,-1); var timesp=ArrayCopy(times).shift();
			for(var i=0;i<timesm.length;i++){
				if(timesm[i]==timesp[i]){
					md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' columns must not contain duplicate timesteps",fieldname)));
					break;
				}
			}
		}
		else{
			md.checkmessage(options.getfieldvalue('message',sprintf("field '%s' should have 2 lines",fieldname)));
		}
	}
}
