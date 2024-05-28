/*
	PointCloudFindNeighborsx.h
*/

#ifndef _POINTCLOUDFLAGNEIGHBORSX_H
#define _POINTCLOUDFLAGNEIGHBORSX_H

#include "../../shared/shared.h"
#include "../../classes/classes.h"

/* local prototypes: */
int PointCloudFindNeighborsx(IssmSeqVec<IssmPDouble>** pflags,double* x, double* y, int nods, double mindistance,double multithread);

/*threading: */
typedef struct{

	double* x;
	double* y;
	int nods;
	double mindistance;
	IssmSeqVec<IssmPDouble>* flags;

} PointCloudFindNeighborsThreadStruct;

void* PointCloudFindNeighborsxt(void* vPointCloudFindNeighborsThreadStruct);

#endif /* _POINTCLOUDFLAGNEIGHBORSX_H */
