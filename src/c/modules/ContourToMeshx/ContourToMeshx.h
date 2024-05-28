/*
	ContourToMeshx.h
*/

#ifndef _CONTOURTOMESHX_H
#define _CONTOURTOMESHX_H

#include "../../shared/shared.h"
#include "../../classes/classes.h"

/*threading: */
typedef struct{

	Contours *contours;
	int       nods;
	int       edgevalue;
	double   *in_nod;
	double   *x;
	double   *y;

} ContourToMeshxThreadStruct;

/* local prototypes: */
int ContourToMeshx(double** pin_nods,double** pin_elem, double* index, double* x, double* y,Contours* contours,char* interptype,int nel,int nods, int edgevalue);

void* ContourToMeshxt(void* vContourToMeshxThreadStruct);

#endif /* _CONTOURTOMESHX_H */
