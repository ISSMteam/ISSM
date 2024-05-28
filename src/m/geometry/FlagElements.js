function FlagElements(md,region){
//FLAGELEMENTS - flag the elements in an region
//
//   The region can be given as a string, or as a javascript array
//
//   Usage: 
//      flag=FlagElements(md,region);
//
//   Example:
//      flag=FlagElements(md,'all');
//      flag=FlagElements(md,'');
//      flag=FlagElements(md,domain);

	//variables
	var flag;
	
	if (typeof region == 'string'){
		if (region === ''){
			flag=NewArrayFill(md.mesh.numberofelements,0);
		}
		else if (region === 'all'){
			flag=NewArrayFill(md.mesh.numberofelements,1);
		}
		else{
			flag=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,region,'element',1);
		}
	}
	else if(IsArray(region)){
		if (region.length==md.mesh.numberofelements){
			flag=region;
		}
		else if (region.length==md.mesh.numberofvertices){
			var flag=NewArrayFill(md.mesh.numberofelements,0);
			for (var i=0;i<md.mesh.numberofelements;i++)
				var sum=0;
				for(var j=0;j<md.mesh.elements[0].length;j++){
					sum += region[md.mesh.element[i][j]-1];
				}
				if (sum==md.mesh.elements[0].length)flag[i]=1;
		}
		else{
			throw Error('Flaglist for region must be of same size as number of elements in model');
		}
	}
	else{
		throw Error('Invalid region option');
	}
	return flag;
}
