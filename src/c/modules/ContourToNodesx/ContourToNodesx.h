/*
	ContourToNodesx.h
*/

#ifndef _CONTOURTONODESX_H
#define _CONTOURTONODESX_H

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

} ContourToNodesThreadStruct;

/* local prototypes: */
int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contour<IssmPDouble>** contours,int numcontours,int edgevalue);
int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contours* contours, int edgevalue);
void* ContourToNodesxt(void* threadstruct);


#endif /* _CONTOURTONODESX_H */
