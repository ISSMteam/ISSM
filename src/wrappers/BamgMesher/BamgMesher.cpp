/*\file BamgMesher.c
 *\brief: mesher that uses the bamg library
 */
#include "./BamgMesher.h"

void BamgMesherUsage(void){/*{{{*/
	_printf0_("\n");
	_printf0_("   usage: [bamgmesh,bamggeom]=" << __FUNCT__ << "(bamgmesh,bamggeom,bamgoptions)\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(BamgMesher_python){

	/*Intermediary*/
	BamgOpts *bamgopts     = NULL;
	BamgMesh *bamgmesh_in  = NULL;
	BamgGeom *bamggeom_in  = NULL;
	BamgMesh *bamgmesh_out = NULL;
	BamgGeom *bamggeom_out = NULL;

	/*Boot module*/
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&BamgMesherUsage);

	/*Initialize outputs*/
	bamggeom_out=new BamgGeom();
	bamgmesh_out=new BamgMesh();

	/*Fetch inputs: */
	FetchData(&bamgopts,BAMGOPTIONS);
	FetchData(&bamggeom_in,BAMGGEOMIN);
	FetchData(&bamgmesh_in,BAMGMESHIN);

	/*Call x layer*/
	Bamgx(bamgmesh_out,bamggeom_out,bamgmesh_in,bamggeom_in,bamgopts);

	/*Generate output Matlab Structures*/
	WriteData(BAMGGEOMOUT,bamggeom_out);
	WriteData(BAMGMESHOUT,bamgmesh_out);

	/*Free resources: */
	delete bamgopts;
	delete bamggeom_in;
	delete bamggeom_out;
	delete bamgmesh_in;
	delete bamgmesh_out;

	/*end module: */
	MODULEEND();
}
