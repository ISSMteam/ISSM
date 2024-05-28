/*!\file BamgConvertMeshx
 */

#include "./BamgConvertMeshx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;

int BamgConvertMeshx(BamgMesh* bamgmesh,BamgGeom* bamggeom,int* index,double* x,double* y,int nods,int nels){

	/*Options*/
	BamgOpts* bamgopts=new BamgOpts();

	/*read mesh*/
	Mesh Th(index,x,y,nods,nels,bamgopts); 

	/*write mesh and geometry*/
	Th.Gh.WriteGeometry(bamggeom,bamgopts);
	Th.WriteMesh(bamgmesh,bamgopts);

	/*clean up and return*/
	delete bamgopts;
	return 1;

}
