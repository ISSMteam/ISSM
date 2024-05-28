function WriteData(fid,prefix){ //{{{
//WRITEDATA - write model field into binary buffer 
//
//   Usage:
//      WriteData(fid,varargin);

	//process options
	var args = Array.prototype.slice.call(arguments);
	var options = new pairoptions(args.slice(2,args.length));
	var name;
	var dataini;
	var data;

	//Get data properties
	if (options.exist('object')){
		//This is a object field, construct string name and data
		obj       = options.getfieldvalue('object');
		fieldname = options.getfieldvalue('fieldname');
		name      = options.getfieldvalue('name',prefix+'.'+fieldname);
		if (options.exist('data')){
			dataini = options.getfieldvalue('data');
		}
		else{
			dataini  = obj[fieldname];
		}
	}
	else{
		//No processing required
		dataini = options.getfieldvalue('data');
		name = options.getfieldvalue('name');
	}
	if (IsArray(dataini)){
	   data=dataini.slice(0);
	}
	else{
		data=dataini;
	}

	format  = options.getfieldvalue('format');
	mattype = options.getfieldvalue('mattype',0);    //only required for matrices
	timeserieslength = options.getfieldvalue('timeserieslength',-1);

	//Scale data if necesarry
	if (options.exist('scale')){
		scale = options.getfieldvalue('scale');
		if (data.length==timeserieslength){
			if (IsArray(data[0])){
				for(var i=0;i<data.length-1;i++){
					for(var j=0;j<data[0].length;j++){
						data[i][j]=scale*data[i][j];
					}
				}
			}
			else{
				for(var i=0;i<data.length-1;i++){
					data[i]=scale*data[i];
				}
			}
		}
		else{
			if (!IsArray(data)) data=data*scale;
			else data=ArrayScale(data,scale);
		}
	}

	if(IsArray(data)){
		if(data.length == timeserieslength){
			var yts = options.getfieldvalue('yts');
			if (IsArray(data[0])){
				for(var j=0;j<data[0].length;j++)data[timeserieslength-1][j]=data[timeserieslength-1][j]*yts;
			}
			else data[timeserieslength-1]=data[timeserieslength-1]*yts;
		}
	}

	let recordlengthtype = 'int';
	if (svnversion >= 22708) {
		recordlengthtype = 'long long';
	}

	//Step 1: write the name to identify this record uniquely
	fid.fwrite(name.length,'int'); 
	fid.fwrite(name,'char'); 

	//Step 2: write the data itself.
	if (format == 'Boolean'){// {{{
		if(IsArray(data)) throw  Error(sprintf("field '%s' cannot be marshalled as it has more than one element!",name));

		//first write length of record
		fid.fwrite(4+4,recordlengthtype);  //1 bool (disguised as an int)+code

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//now write integer
		fid.fwrite(data,'int');  //send an int, not easy to send a bool
	} // }}}
	else if (format == 'Integer'){ // {{{
		if(IsArray(data)) throw  Error(sprintf("field '%s' cannot be marshalled as it has more than one element!",name));

		//first write length of record
		fid.fwrite(4+4,recordlengthtype);  //1 integer + code

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//now write integer
		fid.fwrite(data,'int'); 
	} // }}}
	else if (format == 'Double'){ // {{{
		if(IsArray(data)) throw  Error(sprintf("field '%s' cannot be marshalled as it has more than one element!",name));

		//first write length of record
		fid.fwrite(8+4,recordlengthtype);  //1 double+code

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//now write double
		fid.fwrite(data,'double'); 
	} // }}}
	else if (format == 'String'){ // {{{
		//first write length of record
		fid.fwrite(data.length+4+4,recordlengthtype);  //string + string size + code

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//now write string
		fid.fwrite(data.length,'int'); 
		fid.fwrite(data,'char'); 
	} // }}}
	else if (format == 'BooleanMat'){ // {{{

		//Get size (TODO: use into matlab.js::size)
		var s=[data.length,1]; //vector
		if (IsArray(data[0])) { //matrix
			s[1]=data[0].length;
        } else if (typeof data == 'number') { //scalar
        	s[0]=1; s[1]=1
        } else if (data.length == 0) { //empty matrix/vector
            s[1]=0;
        }

		//if matrix = NaN, then do not write anything
		if (s[0]==1 && s[1]==1 && isNaN(data)){
			s[0]=0; s[1]=0;
		}
		if (typeof data != 'number' && s[0]==1 && s[1]==1 && isNaN(data[0])){
			s[0]=0; s[1]=0;
		}

		//first write length of record
		fid.fwrite(4+4+8*s[0]*s[1]+4+4,recordlengthtype);  //2 integers (32 bits) + the double matrix + code + matrix type

		//write data code and matrix type: 
		fid.fwrite(FormatToCode(format),'int'); 
		fid.fwrite(mattype,'int');

		//now write matrix
		fid.fwrite(s[0],'int'); 
		fid.fwrite(s[1],'int'); 
		if (s[0]*s[1]) fid.fwrite(MatrixToList(data),'double'); //get to the "c" convention, hence the transpose
	} // }}}
	else if (format == 'IntMat'){ // {{{

		//Get size (TODO: use into matlab.js::size)
		var s=[data.length,1]; //vector
		if (IsArray(data[0])) { //matrix
			s[1]=data[0].length;
        } else if (typeof data == 'number') { //scalar
        	s[0]=1; s[1]=1
        } else if (data.length == 0) { //empty matrix/vector
            s[1]=0;
        }

		//if matrix = NaN, then do not write anything
		if (s[0]==1 && s[1]==1 && isNaN(data)){
			s[0]=0; s[1]=0;
		}
		if (typeof data != 'number' && s[0]==1 && s[1]==1 && isNaN(data[0])){
			s[0]=0; s[1]=0;
		}

		//first write length of record
		fid.fwrite(4+4+8*s[0]*s[1]+4+4,recordlengthtype);  //2 integers (32 bits) + the double matrix + code + matrix type

		//write data code and matrix type: 
		fid.fwrite(FormatToCode(format),'int'); 
		fid.fwrite(mattype,'int');

		//now write matrix
		fid.fwrite(s[0],'int'); 
		fid.fwrite(s[1],'int'); 
		if (s[0]*s[1]) fid.fwrite(MatrixToList(data),'double'); //get to the "c" convention, hence the transpose

	} // }}}
	else if (format == 'DoubleMat'){ // {{{

		//Get size (TODO: use into matlab.js::size)
		var s=[data.length,1]; //vector
		if (IsArray(data[0])) { //matrix
			s[1]=data[0].length;
        } else if (typeof data == 'number') { //scalar
        	s[0]=1; s[1]=1
        } else if (data.length == 0) { //empty matrix/vector
            s[1]=0;
        }
        
		//if matrix = NaN, then do not write anything
		if (s[0]==1 && s[1]==1 && isNaN(data)){
			s[0]=0; s[1]=0;
		}
		if (typeof data != 'number' && s[0]==1 && s[1]==1 && isNaN(data[0])){
			s[0]=0; s[1]=0;
		}

		//first write length of record
		var recordlength=4+4+8*s[0]*s[1]+4+4; //2 integers (32 bits) + the double matrix + code + matrix type
		if (recordlength>Math.pow(2,31)) throw Error(sprintf("field '%s' cannot be marshalled because it is larger than 2^31 bytes!",name));
		fid.fwrite(recordlength,recordlengthtype);

		//write data code and matrix type: 
		fid.fwrite(FormatToCode(format),'int'); 
		fid.fwrite(mattype,'int');

		//now write matrix
		fid.fwrite(s[0],'int'); 
		fid.fwrite(s[1],'int'); 
		if (s[0]*s[1]) fid.fwrite(MatrixToList(data),'double'); //get to the "c" convention, hence the transpose
	} // }}}
	else if (format == 'MatArray'){ // {{{

		numrecords=data.length;

		//first get length of record
		recordlength=4+4; //number of records + code
		for (var i=0;i<numrecords;i++){
			matrix=data[i];
			var s=[matrix.length,1];
			if(IsArray(matrix[0]))s[1]=matrix[0].length;

			recordlength=recordlength+4*2+ //row and col of matrix
				s[0]*s[1]*8; //matrix of doubles
		}

		//write length of record
		fid.fwrite(recordlength,recordlengthtype); 

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//write data, first number of records
		fid.fwrite(numrecords,'int'); 

		//write each matrix: 
		for (var i=0;i<numrecords;i++){
			matrix=data[i];
			var s=[matrix.length,1];
			if(IsArray(matrix[0]))s[1]=matrix[0].length;

			fid.fwrite(s[0],'int'); 
			fid.fwrite(s[1],'int'); 
			fid.fwrite(MatrixToList(matrix),'double');
		}
	} // }}}
	else if (format == 'StringArray'){ // {{{

		//first get length of string array: 
		num=data.length;
		if ((typeof data[0] == 'number') & num==1 & isNaN(data[0])){
			num = 0;
		}

		//now get length of record: 
		recordlength=4+4; //for length of array + code
		for (var i=0;i<num;i++){
			string=data[i];
			recordlength=recordlength+4+string.length; //for each string
		}

		//write length of record
		fid.fwrite(recordlength,recordlengthtype); 

		//write data code: 
		fid.fwrite(FormatToCode(format),'int'); 

		//now write length of string array
		fid.fwrite(num,'int'); 

		//now write the strings
		for (var i=0;i<num;i++){
			string=data[i];
			fid.fwrite(string.length,'int'); 
			fid.fwrite(string,'char'); 
		}
	} // }}}
	else { 
		throw Error(sprintf("WriteData error message: data type: %s not supported yet! ('%s')",
					format.toString(),name));
	}
} //}}}
function FormatToCode(format){ // {{{
	//This routine takes the format string, and hardcodes it into an integer, which 
	//is passed along the record, in order to identify the nature of the dataset being 
	//sent.
	if  (format == 'Boolean') code=1;
	else if (format == 'Integer') code=2;
	else if (format == 'Double') code=3;
	else if (format == 'String') code=4;
	else if (format == 'BooleanMat') code=5;
	else if (format == 'IntMat') code=6;
	else if (format == 'DoubleMat') code=7;
	else if (format == 'MatArray') code=8;
	else if (format == 'StringArray') code=9;
	else throw Error('FormatToCode error message: data type not supported yet!');
	return code;
}// }}}
