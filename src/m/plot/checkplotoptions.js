function checkplotoptions(md,options){ //{{{
	//PARSE_OPTIONS - build a structure that holds all plot options
	//
	//   Usage:
	//      checkplotoptions(md,options);
	//
	//   See also: PLOTMODEL

	//units
	if (options.exist('unit')){
		if (options.getfieldvalue('unit') === 'km'){
			options.changefieldvalue('unit',Math.pow(10,-3));
		}
		if (options.getfieldvalue('unit') === '100km'){
			options.changefieldvalue('unit',Math.pow(10,-5));
		}
	}

	//density
	if (options.exist('density')){
		density=options.getfieldvalue('density');
		options.changefieldvalue('density',Math.abs(Math.ceil(density)));
	}

	//Show section
	if (options.exist('showsection')){
		if (options.getfieldvalue('showsection') === 'on'){
			options.changefieldvalue('showsection',4);
		}
	}

	//smooth values
	if (options.exist('smooth')){
		if (options.getfieldvalue('smooth') === 'on'){
			options.changefieldvalue('smooth',0);
		}
	}

	//contouronly values
	if (options.exist('contouronly')){
		if (options.getfieldvalue('contouronly') === 'on'){
			options.changefieldvalue('contouronly',1);
		}
	}

	//Colorbar;
	if (options.exist('colorbar')){
		if (options.getfieldvalue('colorbar') === 'on'){
			options.changefieldvalue('colorbar',1);
		}
		else if (options.getfieldvalue('colorbar') === 'off'){
			options.changefieldvalue('colorbar',0);
		}
	}

	//text
	if (options.exist('text')){

		//1: textvalue
		textvalues=options.getfieldvalue('text');

		//ischar if only one expstyle -> create a cell
		if (typeof textvalues === 'string'){
			textvalues=[textvalues];
			numtext=1;
		}
		else if (IsArray(textvalues)){
			numtext=textvalues.length;
		}
		else throw Error("plot error message: ''text'' option should be either a string or a cell");

		//2: textweight
		if (options.exist('textweight')){

			textweightvalues=options.getfieldvalue('textweight');

			//ischar if only one textweight -> create a cell
			if (typeof textweightvalues === 'string'){
				textweightvalues=[textweightvalues];
			}
			else if (!IsArray(textweightvalues)){
				throw Error("plot error message: ''textweight'' option should be either a string or a cell");
			}
		}
		else{
			textweightvalues=['n'];
		}
		if (textweightvalues.length==1){
			var value=textweightvalues[0];
			for (var i=0;i<numtext-1;i++)textweightvalues.push(value);
		}

		//3: textsize
		if (options.exist('textsize')){
			textsizevalues=options.getfieldvalue('textsize');
		}
		//ischar if only one textsize -> create a cell
		if (typeof textsizevalues === 'number'){
			textsizevalues=[textsizevalues];
		}
		else if (!IsArray(textsizevalues)){
			throw Error("plot error message: ''textsize'' option should be either a number or a cell");
		}
		else{
			textsizevalues=[14];
		}
		if (textsizevalues.length==1){
			var value=textsizevalues[0];
			for (var i=0;i<numtext-1;i++)textsizevalues.push(value);
		}
			
		//4: textcolor
		if (options.exist('textcolor')){
			textcolorvalues=options.getfieldvalue('textcolor');
		}
		if (typeof textcolorvalues === 'string'){ //ischar if only one textcolor -> create a cell
			textcolorvalues=[textcolorvalues];
		}
		else if (!IsArray(textcolorvalues)){
			throw Error("plot error message: ''textcolor'' option should be either a string or a cell");
		}
		else textcolorvalues=['k'];

		if (textcolorvalues.length==1){
			var value=textcolorvalues[0];
			for (var i=0;i<numtext-1;i++)textcolorvalues.push(value);
		}
		
		//5: textposition
		if (options.exist('textposition')){
			textpositionvalues=options.getfieldvalue('textposition');
		}
		//ischar if only one textposition -> create a cell
		if (typeof textpositionvalues === 'number'){
			textpositionvalues=[textpositionvalues];
		}
		else if(!IsArray(textpositionvalues)){
			throw Error("plot error message: ''textposition'' option should be either a string or a cell");
		}
		else throw Error("plot error message: ''textposition'' option is missing");
			
		//6: textrotation
		if (options.exist('textrotation')){
			textrotationvalues=options.getfieldvalue('textrotation');
		}
		//ischar if only one textsize -> create a cell
		if (typeof textrotationvalues === 'number'){
			textrotationvalues=[textrotationvalues];
		}
		else if (!IsArray(textrotationvalues)){
			throw Error("plot error message: ''textrotation'' option should be either a number or a cell");
		}
		else textrotationvalues=[0];
		
		if (textrotationvalues.length==1){
			var value=textrotationvalues[0];
			for (var i=0;i<numtext-1;i++)textrotationvalues.push(value);
		}
			
		options.changefieldvalue('text',textvalues);
		options.changefieldvalue('textsize',textsizevalues);
		options.changefieldvalue('textweight',textweightvalues);
		options.changefieldvalue('textcolor',textcolorvalues);
		options.changefieldvalue('textposition',textpositionvalues);
		options.changefieldvalue('textrotation',textrotationvalues);
	}

	//expdisp
	expdispvaluesarray=[];
	expstylevaluesarray=[];
	expstylevalues=[];
	if (options.exist('expstyle')){
		expstylevalues=options.getfieldvalue('expstyle');
		//ischar if only one expstyle -> create a cell
		if (typeof expstylevalues === 'string'){
			expstylevalues=[expstylevalues];
		}
		options.changefieldvalue('expdisp',expdispvaluesarray);
	}
		
	if (options.exist('expdisp')){
		expdispvalues=options.getfieldvalue('expdisp');
	
		//ischar if only one expstyle -> create a cell
		if (typeof expdispvalues === 'string'){
			expdispvalues=[expdispvalues];
		}
		for (var i=0; i< expdispvalues.length;i++){
			expdispvaluesarray.push(expdispvalues[i]);
			if (expstylevalues.length>i){
				expstylevaluesarray.push(expstylevalues[i]);
			}
			else{
				expstylevaluesarray.push('g-');
			}
		}
		options.changefieldvalue('expstyle',expstylevaluesarray);
	}

	//latlonnumbering
	if (options.exist('latlonclick')){
		if (options.getfieldvalue('latlonclick') === 'on'){
			options.changefieldvalue('latlonclick',1);
		}
	}

	//north arrow
	if (options.exist('northarrow')){
	   if (options.getfieldvalue('northarrow') === 'on'){
		   
		   //default values
		   Lx=ArrayMax(md.mesh.y)-ArrayMin(md.mesh.y);
		   Ly=ArrayMax(md.mesh.y)-ArrayMin(md.mesh.y);
		  
		   //default values
		   options.changefieldvalue('northarrow',[ArrayMin(md.mesh.x)+1/6*Lx,ArrayMin(md.mesh.y)+5/6*Ly,1/15*Ly,0.25,1/250*Ly]);
	   }
	}

	//scale ruler
	if (options.exist('scaleruler')){
	   if (options.getfieldvalue('scaleruler') === 'on'){
		   //default values
		   Lx=ArrayMax(md.mesh.x)-ArrayMin(md.mesh.x);
		   Ly=ArrayMax(md.mesh.y)-ArrayMin(md.mesh.y);
		   
		   //default values
		   options.changefieldvalue('scaleruler',[ArrayMin(md.mesh.x)+6/8*Lx, ArrayMin(md.mesh.y)+1/10*Ly, Math.pow(10,(Math.ceil(Math.log10(Lx))))/5, Math.floor(Lx/100), 5]);
	   }
	}
} //}}}
