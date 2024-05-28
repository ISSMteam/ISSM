/*\file BamgConvertMesh.c
 *\brief: bamg module.
 */
#include "./BamgConvertMesh.h"

void BamgConvertMeshUsage(void){/*{{{*/
	_printf0_("BAMGCONVERTMESH - convert [x y index] to a bamg geom and mesh geom");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      [bamggeom bamgmesh]=BamgConvertMesh(index,x,y)\n");
	_printf0_("      index: index of the mesh\n");
	_printf0_("      x,y: coordinates of the nodes\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(BamgConvertMesh_python){

	/*input: */
	int    *index      = NULL;
	double *x          = NULL;
	double *y          = NULL;
	int     nods,nels,test1,test2;

	/*Output*/
	BamgMesh *bamgmesh = NULL;
	BamgGeom *bamggeom = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&BamgConvertMeshUsage);

	/*Initialize Bamg outputs*/
	bamggeom=new BamgGeom();
	bamgmesh=new BamgMesh();

	/*Input datasets: */
	FetchData(&index,&nels,&test1,INDEXHANDLE);
	FetchData(&x,&nods,XHANDLE);
	FetchData(&y,&test2,YHANDLE);

	/*Check inputs*/
	if(nels<0) _error_("Number of elements must be positive, check index number of lines");
	if(nods<0) _error_("Number of nods must be positive, check x and y sizes");
	if(test1!=3) _error_("index should have 3 columns");
	if(test2!=nods) _error_("x and y do not have the same length");

	/* Run core computations: */
	BamgConvertMeshx(bamgmesh,bamggeom,index,x,y,nods,nels);

	/*Generate output Matlab Structures*/
	WriteData(BAMGGEOMOUT,bamggeom);
	WriteData(BAMGMESHOUT,bamgmesh);

	/*Clean up*/
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	delete bamggeom;
	delete bamgmesh;

	/*end module: */
	MODULEEND();
}
