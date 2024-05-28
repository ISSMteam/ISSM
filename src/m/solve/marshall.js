function marshall(md){
//MARSHALL - outputs a typed array buffer to be send to the issm module.
//
//   The routine creates a compatible binary stream from @model md
//   This binary stream will be used for single cpu runs using the issm module.
//
//   Usage:
//      fid=marshall(md)

	if (md.verbose.solution){
		console.log('marshalling file ' + md.miscellaneous.name + '.bin');
	}

	//open file for binary writing
	var fid=new fileptr('mode','w');

	//Go through all model fields: check that it is a class and call checkconsistency
	for (field in md){

		//Some properties do not need to be marshalled
		if (field == 'results' | field =='radaroverlay' | field == 'toolkits' | field =='cluster' | field == 'priv') continue;
		
		//Check that current field is a class
		if(typeof md[field] == 'function'){
			continue;
		}

		//Marshall current object
		md[field].marshall(md,['md.'+field],fid);
	}

	//Last, write "md.EOF" to make sure that the binary file is not corrupt
	WriteData(fid,'XXX','name','md.EOF','data',true,'format','Boolean');
	return fid;
}
