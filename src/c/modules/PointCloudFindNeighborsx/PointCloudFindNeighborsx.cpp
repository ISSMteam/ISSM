/*! \file  PointCloudFindNeighborsx.c
 */

#include "./PointCloudFindNeighborsx.h"

int PointCloudFindNeighborsx(IssmSeqVec<IssmPDouble>** pflags,double* x, double* y, int nods, double mindistance,double multithread){

	/*output: */
	IssmSeqVec<IssmPDouble>* flags=NULL;
	flags=new IssmSeqVec<IssmPDouble>(nods);

	/*threading: */
	int num=_NUMTHREADS_;
	if(!multithread)num=1;

	/*initialize thread parameters: */
	PointCloudFindNeighborsThreadStruct gate;
	gate.x           = x;
	gate.y           = y;
	gate.nods        = nods;
	gate.mindistance = mindistance;
	gate.flags       = flags;

	/*launch the thread manager with InterpFromGridToMeshxt as a core: */
	LaunchThread(PointCloudFindNeighborsxt,(void*)&gate,num);

	/*Assemble vector: */
	flags->Assemble();

	/*Assign output pointers: */
	*pflags=flags;

	return 1;
}
