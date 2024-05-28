function SetIceShelfBC(md) { 
//SETICESHELFBC - Create the boundary conditions for stressbalance and thermal models for a  Ice Shelf with Ice Front
//
//   Neumann BC are used on the ice front (an ANRGUS contour around the ice front
//   must be given in input)
//   Dirichlet BC are used elsewhere for stressbalance
//
//   Usage:
//      md=SetIceShelfBC(md,varargin)
//
//   Example:
//      SetIceShelfBC(md);
//      SetIceShelfBC(md,'Front.exp');
//
//   See also: SETICESHEETBC, SETMARINEICESHEETBC

	//node on Dirichlet (boundary and ~icefront)
	if (arguments.length==2){
		icefront=arguments[1];
		nodeinsideicefront=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,icefront,'node',2);
		nodeonicefront=ArrayAnd(md.mesh.vertexonboundary,nodeinsideicefront);
	}
	else if(arguments.length==1){
		nodeonicefront=NewArrayFill(md.mesh.numberofvertices,0);
	}
	else{
		throw Error('SetIceShelfBC usage error');
	}

	md.stressbalance.spcvx=NewArrayFill(md.mesh.numberofvertices,NaN); 
	md.stressbalance.spcvy=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.stressbalance.spcvz=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.stressbalance.referential=Create2DArray(md.mesh.numberofvertices,6);
	for(var i=0;i<md.mesh.numberofvertices;i++)for(var j=0;j<6;j++)md.stressbalance.referential[i][j]=NaN;
	md.stressbalance.loadingforce=Create2DArray(md.mesh.numberofvertices,3);
	for(var i=0;i<md.mesh.numberofvertices;i++)for(var j=0;j<3;j++)md.stressbalance.loadingforce[i][j]=0;

	//Ice front position: 
	pos=ArrayFind(nodeonicefront,1);
	for(var i=0;i<pos.length;i++)md.mask.ice_levelset[pos[i]]=0;

	//First find segments that are not completely on the front
	if (md.mesh.elementtype() === 'Penta'){
		numbernodesfront=4;
	}
	else if (md.mesh.elementtype() === 'Tria'){
		numbernodesfront=2;
	}
	else{
		throw Error('mesh type not supported yet');
	}
	var obs=false;
	if((md.inversion.vx_obs.length == md.mesh.numberofvertices) & (md.inversion.vy_obs.length==md.mesh.numberofvertices))obs=true;

	if(obs==true){
		console.log('      boundary conditions for stressbalance model: setting spc as observed velocities');
	}
	else{
		console.log('      boundary conditions for stressbalance model: setting spc as zero');
	}
	for(var i=0;i<md.mesh.segments.length;i++){
		var sum=0;
		for (var j=0;j<numbernodesfront;j++) sum+=md.mask.ice_levelset[md.mesh.segments[i][j]-1];
		if(sum!=0){
			for (var j=0;j<numbernodesfront;j++){
				if(obs==false){
					md.stressbalance.spcvx[md.mesh.segments[i][j]-1]=0;
					md.stressbalance.spcvy[md.mesh.segments[i][j]-1]=0;
				}
				else{
					md.stressbalance.spcvx[md.mesh.segments[i][j]-1]=md.inversion.vx_obs[md.mesh.segments[i][j]-1];
					md.stressbalance.spcvy[md.mesh.segments[i][j]-1]=md.inversion.vy_obs[md.mesh.segments[i][j]-1];
				}
				md.stressbalance.spcvz[md.mesh.segments[i][j]-1]=0;

			}
		}
	}

	//Initialize surface and basal forcings
	md.smb.initialize(md);
	md.basalforcings.initialize(md);

	//Deal with other boundary conditions
	if (isNaN(md.balancethickness.thickening_rate)){
		md.balancethickness.thickening_rate=NewArrayFill(md.mesh.numberofvertices,0);
		console.log('      no balancethickness.thickening_rate specified: values set as zero');
	}
		
	md.masstransport.spcthickness=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.balancethickness.spcthickness=NewArrayFill(md.mesh.numberofvertices,NaN);
	md.damage.spcdamage=NewArrayFill(md.mesh.numberofvertices,NaN);

	if (md.initialization.temperature.length==md.mesh.numberofvertices){
		md.thermal.spctemperature=NewArrayFill(md.mesh.numberofvertices,NaN);
		if ('vertexonsurface' in md.mesh){
			pos=ArrayFind(md.mesh.vertexonsurface,1);
			for(var i=0;i<pos.length;i++)md.thermal.spctemperature[i]=md.initialization.temperature[i]; //impose observed temperature on surface
		}
		if (md.basalforcings.geothermalflux.length != md.mesh.numberofvertices){
			md.basalforcings.geothermalflux=NewArrayFill(md.mesh.numberofvertices,0);
		}
	}
	else{
		console.log('      no thermal boundary conditions created: no observed temperature found');
	}
}
