function setmask(md,floatingice,groundedice){
//SETMASK - establish boundaries between grounded and floating ice.
//
//   By default, ice is considered grounded. The contour floatingice defines nodes 
//   for which ice is floating. The contour groundedice defines nodes inside a floatingice, 
//   that are grounded (ie: ice rises, islands, etc ...)
//   All inputs are either strings or actually javascript arrays (included in the html file)
//   For example: 
//
//	   floatingice[0]['x']=[0,0,0,1];
//	   floatingice[0]['y']=[0,1,1,1];
//	   floatingice[1]['x']=[0,0.5,0.5,.5];
//	   floatingice[1]['y']=[0,.5,.5,.5];
//
//
//   Usage:
//      md=setmask(md,floatingice,groundedice)
//
//   Examples:
//      md=setmask(md,'all','');
//      md=setmask(md,iceshelves,islands);

	//variables: 
	var  icedomain=[];
	
	//some checks on list of arguments
	if (!((arguments.length==3) | (arguments.length==5))){
		throw Error('mask error message: wrong usage.');
	}

	if(arguments.length>3){
		if (arguments[3]=='icedomain'){
			icedomain=arguments[4];
		}
		else{
			throw Error('mask error message: wrong field specified. Only icedomain allowed for now.');
		}
		if (IsArray(icedomain)){
			throw Error('setmask error message: icedomain should be an array!');
		}
	}
	
	//Get assigned fields
	var x=md.mesh.x;
	var y=md.mesh.y;
	var elements=md.mesh.elements;

	//Assign elementonfloatingice, elementongroundedice, vertexongroundedice and vertexonfloatingice. 
	//Only change at your own peril! This is synchronized heavily with the GroundingLineMigration module. 
	elementonfloatingice=FlagElements(md,floatingice);
	elementongroundedice=FlagElements(md,groundedice);

	//Because groundedice nodes and elements can be included into an floatingice, we need to update. Remember, all the previous 
	//arrays come from domain outlines that can intersect one another: 
	elementonfloatingice=ArrayAnd(elementonfloatingice,ArrayNot(elementongroundedice));
	elementongroundedice=ArrayNot(elementonfloatingice);

	//the order here is important. we choose vertexongroundedice as default on the grounding line.
	vertexonfloatingice=NewArrayFill(md.mesh.numberofvertices,0);
	vertexongroundedice=NewArrayFill(md.mesh.numberofvertices,0);
	pos=ArrayFind(elementongroundedice,1); for (var i=0;i<pos.length;i++)for(var j=0;j<3;j++) vertexongroundedice[md.mesh.elements[i,j]-1]=1;
	pos=ArrayFind(vertexongroundedice,0); for (var i=0;i<pos.length;i++)vertexonfloatingice[i]=1;

	//level sets
	ocean_levelset=vertexongroundedice;
	pos=ArrayFind(vertexongroundedice,0);for(var i=0;i<pos.length;i++) ocean_levelset[i]=-1;
	md.mask.ocean_levelset=ocean_levelset;

	if(arguments.length>3){
		md.mask.ice_levelset = NewArrayFill(md.mesh.numberofvertices,1.0);
		//use contourtomesh to set ice values inside ice domain
		//[vertexinsideicedomain,elementinsideicedomain]=ContourToMesh(elements,x,y,icedomain,'node',1);
		pos=ArrayFind(vertexinsideicedomain,1.0);for(var i=0;i<pos.length;i++) md.mask.ice_levelset[pos]=-1;
	}
	else{
		md.mask.ice_levelset = NewArrayFill(md.mesh.numberofvertices,-1);
	}
}
