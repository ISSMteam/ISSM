function processdata(md,data,options){ //{{{
	//PROCESSDATA - process data to be plotted
	//
	//   datatype = 1 -> elements
	//   datatype = 2 -> nodes
	//   datatype = 3 -> node quivers
	//   datatype = 4 -> patch
	//   datatype = 5 -> nodes transient
	//
	//   Usage:
	//      var array = processdata(md,data,options);
	//      var data = array[0]; 
	//      var datatype = array[1];
	//
	//   See also: PLOTMODEL, PROCESSMESH

	//check format
	if ( data.length ==0 | data === [] | typeof data === 'number' | ArrayAnyNaN(data) ){
		throw Error('plotmodel error message: data provided is empty');
	}

	//Process NaN if any (do not know before mask is applied)
	if (options.exist('nan')){
		var valuefornan=options.getfieldvalue('nan',0);
		for (var i=0;i<data.length;i++)if(IsNaN(data[i]))data[i]=valuefornan;
	}

	//special case for mesh 2dvertical
	if (md.mesh.domaintype() === '2Dvertical'){
		return processdata(md.mesh,md,data,options);
	}

	//needed later on
	var numberofvertices2d, numberofelements2d;
	if ('numberofvertices2d' in md.mesh){
		numberofvertices2d=md.mesh.numberofvertices2d; 
		numberofelements2d=md.mesh.numberofelements2d; 
	}
	else {
		numberofvertices2d=NaN;
		numberofelements2d=NaN;
	}

	//initialize datatype
	var datatype=0;

	//get datasize
	var datasize=data.length;

	//transpose data if necessary
	if (data[0].length > data.length){
		throw Error('processdata error message: you need to tranpose your data!');
	}


	//check length
	if (datasize != md.mesh.numberofvertices & datasize !=md.mesh.numberofelements & datasize!=md.mesh.numberofvertices*6 & 
			((md.mesh.domaintype() === '3D') & !(datasize==numberofelements2d | datasize==numberofvertices2d))){
		throw Error('plotmodel error message: data not supported yet');
	}


	//quiver?
	if (Array.isArray(data[0])){
		datatype=3;

		//check number of columns, add zeros if necessary,
		if (md.mesh.dimension()==3){
			if (data[0].length==2){
				data=[data, NewArrayFill(data.length,1)];
			}
			else if (data[0].length!=3){
				throw Error('plotmodel error message: data provided should have 2 or 3 columns for quiver plot, and 1 for regular plot');
			}
		}
	}

	//treat the case datasize(1)=6*nodes
	if (datasize==6*md.mesh.numberofvertices){
		//keep the only norm of data
		data1=new Array(md.mesh.numberofvertices);
		data2=new Array(md.mesh.numberofvertices);
		data=new Array(md.mesh.numberofvertices);
		for(var i=0;i<md.mesh.numberofvertices;i++){
			data1[i]=data[6*i+0];
			data2[i]=data[6*i+1];
			data[i]=Math.sqrt(pow(data1[i],2),pow(data2[i],2));
		}
		datasize=md.mesh.numberofvertices;
		//---> go to node data
	}

	//treat the case datasize(1)=nodes2d
	if (md.mesh.dimension()==3 & datasize==numberofvertices2d){
		data=project3d(md,'vector',data,'type','node');
		datasize=md.mesh.numberofvertices;
		//---> go to node data
	}

	//treat the case datasize=nodes2d
	if (md.mesh.dimension()==3 & datasize==numberofelements2d){
		data=project3d(md,'vector',data,'type','element');
		datasize=md.mesh.numberofelements;
		//---> go to node data
	}

	//smoothing?
	if (options.exist('smooth')){
		data=averaging(md,data,options.getfieldvalue('smooth'));
		datasize=md.mesh.numberofvertices;
		//---> go to node data
	}

	//element data
	if (datasize==md.mesh.numberofelements & !Array.isArray(data[0])){

		//Initialize datatype if non patch
		if(datatype!=4 & datatype!=5){
			datatype=1;
		}

		//Mask?
		if(options.exist('mask')){
			flags=options.getfieldvalue('mask');
			if(flags.length==md.mesh.numberofvertices){
				for(var i=0;i<md.mesh.numberofelements;i++){
					var nanify=0;
					for(var j=0;j<md.mesh.elements[0].length;j++){
						if (flags[md.mesh.elements[i][j]-1]==0)nanify=1;
					}
					if(nanify) for(var j=0;j<md.mesh.elements[0].length;j++)data[md.mesh.elements[i][j]-1]=NaN;
				}
			}
			else if (flags.length==md.mesh.numberofelements){
				for(var i=0;i<md.mesh.numberofelements;i++)if (flags[i]==0)data[i]=NaN;
			}
			else{
				console.log('plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements)');
			}
		}

		//log?
		if (options.getfieldvalue('log','off')!='off'){
			var bounds=options.getfieldvalue('caxis',[ArrayMin(data),ArrayMax(data)]);
			for(var i=0;i<md.mesh.numberofelements;i++)if(data[i]<bounds[0])data[i]=bounds[0];
			for(var i=0;i<md.mesh.numberofelements;i++)if(data[i]<=0){
				throw Error("Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])");
			}
			for(var i=0;i<md.mesh.numberofelements;i++){
				if(!IsNaN(data[i])){
					data[i]=Math.log10(data[i])/Math.log10(options.getfieldvalue('log',10));
				}
			}
		}
	}

	//node data
	if (datasize==md.mesh.numberofvertices & !Array.isArray(data[0])){
		datatype=2;

		//Mask?
		if (options.exist('mask')){
			flags=options.getfieldvalue('mask');
			if (flags.length==md.mesh.numberofvertices){
				for(var i=0;i<md.mesh.numberofvertices;i++){
					if(flags[i]==0)data[i]=NaN;
				}
			}
			else if( length(flags)==md.mesh.numberofelements){
				for(var i=0;i<md.mesh.numberofelements;i++){
					if(flags[i]==0){
						for(var j=0;j<md.mesh.elements[0].length;j++){
							data[md.mesh.elements[i][j]-1]=NaN;
						}
					}
				}
			}
			else{
				console.log("plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements");
			}
		}

		//log?
		if (options.getfieldvalue('log','off')!='off'){
			var bounds=options.getfieldvalue('caxis',[ArrayMin(data),ArrayMax(data)]);
			for(var i=0;i<md.mesh.numberofvertices;i++)if(data[i]<bounds[0])data[i]=bounds[0];
			for(var i=0;i<md.mesh.numberofvertices;i++)if(data[i]>bounds[1])data[i]=bounds[1];
			for(var i=0;i<md.mesh.numberofvertices;i++)if(data[i]<=0){
				throw Error("Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])");
			}
			for(var i=0;i<md.mesh.numberofvertices;i++){
			   data[i]=Math.log10(data[i])/Math.log10(options.getfieldvalue('log',10));
			}
		}
	}
	
	//node transient data
    if (datasize==md.mesh.numberofvertices+1){
        datatype=5;
		
		//log?	
		if (options.getfieldvalue('log','off')!='off'){
			var bounds=options.getfieldvalue('caxis',[ArrayMin(data),ArrayMax(data)]);
			for(var i=0;i<md.mesh.numberofvertices;i++) {
				for(var j=0;j<data[i].length;j++) {
					if(data[i][j]<bounds[0])data[i][j]=bounds[0];
				}
			}
			for(var i=0;i<md.mesh.numberofvertices;i++) {
				for(var j=0;j<data[i].length;j++) {
					if(data[i][j]>bounds[1])data[i][j]=bounds[1];
				}
			}
			for(var i=0;i<md.mesh.numberofvertices;i++) {
				for(var j=0;j<data[i].length;j++) {
					if(data[i][j]<=0) {
						throw Error("Log option cannot be applied on negative values. Use caxis option (Rignot''s settings: [1.5 max(data)])");
					}
				}
			}
			for(var i=0;i<md.mesh.numberofvertices;i++){
				for(var j=0;j<data[i].length;j++) {
					data[i][j]=Math.log10(data[i][j])/Math.log10(options.getfieldvalue('log',10));
				}
			}
		}
    }
	

	//layer projection? 
	if (options.getfieldvalue('layer',0)>=1){
		data=project2d(md,data,options.getfieldvalue('layer')); //project onto 2d mesh
	}

	//control arrow density if quiverplot: not done yet since conversion of matlab to javascript.
	/*if (datatype==3 & options.exist('density')){
		databak=data;
		data=NewArrayFill(datasize,NaN);

		density=options.getfieldvalue('density');
		data(1:density:end,:)=databak(1:density:end,:);
		clear databak
	}*/

	/*if (datatype==3){ //not done yet
		//Mask?
		if (options.exist('mask')){
			flags=options.getfieldvalue('mask');
			pos=find(~flags);
			if (flags.length==md.mesh.numberofvertices){
			   data(pos,:)=NaN;
			}
			else if (flags.length==md.mesh.numberofelements){
				data(md.mesh.elements(pos,:),:)=NaN;
			}
			else{
				console.log("plotmodel warning: mask length not supported yet (supported length are md.mesh.numberofvertices and md.mesh.numberofelements");
			}
		}
	}*/

	//OK, if datatype=0 error out
	if (datatype==0){
	   throw Error('data provided not recognized or not supported');
	}

	return [data,datatype];
} //}}}
