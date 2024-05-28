function parseresultsfrombuffer(md,buffer,buffersize){ //{{{

	//Open file
	results=[];
	
	var fid = new fileptr('mode','r','buffer',buffer,'buffersize',buffersize);

	//Read fields until the end of the file.
	result  = ReadData(fid,md);

	if (MapIsEmpty(result))throw Error('no results found in binary memory buffer');
	
	var counter = 0;
	var check_nomoresteps=0;
	var step    = result['step'];
	var index;

	while (!MapIsEmpty(result)){

		if (check_nomoresteps){
			//check that the new result does not add a step, which would be an error: 
			if (result['step']>=1)throw Error('parsing results for a steady-state core, which incorporates transient results!');
		}

		//Check step, increase counter if this is a new step
		if(step!=result['step'] & result['step']>1){
			counter = counter + 1;
			step    = result['step'];
		}

		//Add result
		if(result['step']==0){
			//if we have a step = 0, this is a steady state solutoin, don't expect more steps. 
			index = 0;
			check_nomoresteps=1;
		}
		else if(result['step']==1){
			index = 0;
		}
		else index = counter;

		if(index>results.length-1)results.push({});
		for(var i=results.length-1;i<index-1;i++)results[i]={};
		results[index][result['fieldname']]=result['field'];
		
		//Get time and step
		if(result['time']!=-9999){
			results[index]['step']=result['step'];
		}
		if(result['time']!=-9999){
			results[index]['time']=result['time'];
		}

		//read next result
		result  = ReadData(fid,md);
	}
	return results;
} // }}}
function ReadData(fid,md){ //{{{
//READDATA - ...
//
//   Usage:
//      field=ReadData(fid,md)

	//read field
	var length=fid.fread(1,'int');
	
	var result={};

	if (length!==-1){
		fieldname=fid.fread(length,'char');
		time=fid.fread(1,'double');
		step=fid.fread(1,'int');
		type=fid.fread(1,'int');
		M=fid.fread(1,'int');

		if (type==1) field=fid.fread(M,'double');
		else if (type==2) field=fid.fread(M,'char');
		else if (type==3) {
			N=fid.fread(1,'int');
			field=fid.fread(N*M,'double');
		}
		else if (type==4) {
			N=fid.fread(1,'int');
			field=fid.fread(N*M,'int');
		}
		else throw Error(sprintf("%s%i",'ReadData error message: cannot read data of type ',type));

		//Process units here FIXME: this should not be done here!
		var yts=md.constants.yts;
		if (fieldname == 'BalancethicknessThickeningRate') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'HydrologyWaterVx') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'HydrologyWaterVy') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Vx') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Vy') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Vz') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Vel') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'BasalforcingsGroundediceMeltingRate') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'BasalforcingsFloatingiceMeltingRate') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'TotalSmb') for (var i=0;i<field.length;i++)field[i]= field[i]/Math.pow(10,12)*yts; //(GigaTon/year)
		else if (fieldname == 'TotalSmbScaled') for (var i=0;i<field.length;i++)field[i]= field[i]/Math.pow(10,12)*yts; //(GigaTon/year)
		else if (fieldname == 'GroundinglineMassFlux') for (var i=0;i<field.length;i++)field[i]= field[i]/Math.pow(10,12)*yts; //(GigaTon/year)
		else if (fieldname == 'IcefrontMassFlux') for (var i=0;i<field.length;i++)field[i]= field[i]/Math.pow(10,12)*yts; //(GigaTon/year)
		else if (fieldname == 'IcefrontMassFluxLevelset') for (var i=0;i<field.length;i++)field[i]= field[i]/Math.pow(10,12)*yts; //(GigaTon/year)
		else if (fieldname == 'SmbMassBalance') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'SmbPrecipitation') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'SmbRunoff') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'SmbCondensation') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'SmbAccumulation') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'SmbMelt') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'CalvingCalvingrate') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Calvingratey') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'Calvingratex') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;
		else if (fieldname == 'CalvingMeltingrate') for (var i=0;i<field.length;i++)field[i]= field[i]*yts;

		result['fieldname']=fieldname;
		result['time']=time;
		if (result['time']!=-9999) result['time']=time/yts;
		result['step']=step;
		result['field']=field;
	}
	return result;

} // }}}
