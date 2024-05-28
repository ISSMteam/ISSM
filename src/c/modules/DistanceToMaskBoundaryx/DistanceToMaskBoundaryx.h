/*
	DistanceToMaskBoundaryx.h
*/

#ifndef _DISTANCETOMASKBOUNDARYX_H
#define _DISTANCETOMASKBOUNDARYX_H

#include "../../shared/shared.h"
#include "../../classes/classes.h"

/*threading: */
typedef struct{

	int       nods;
	IssmDouble   *distance;
	IssmDouble   *x;
	IssmDouble   *y;
	IssmDouble   *mask;

} DistanceToMaskBoundaryxThreadStruct;

/* local prototypes: */
int DistanceToMaskBoundaryx(IssmDouble** pdistance,IssmDouble* x, IssmDouble* y, IssmDouble* mask, int nods);

void* DistanceToMaskBoundaryxt(void* vDistanceToMaskBoundaryxThreadStruct);

#endif /* _DISTANCETOMASKBOUNDARYX_H */
