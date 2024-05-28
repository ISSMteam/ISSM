/*
	ExpToLevelSetx.h
*/

#ifndef _EXPTOLEVELSETX_H
#define _EXPTOLEVELSETX_H

#include "../../shared/shared.h"
#include "../../classes/classes.h"

/*threading: */
typedef struct{

	Contours *contours;
	int       nods;
	double   *distance;
	double   *x;
	double   *y;

} ExpToLevelSetxThreadStruct;

/* local prototypes: */
int ExpToLevelSetx(double** pdistance,double* x, double* y,int nods, Contours* contours);

void* ExpToLevelSetxt(void* vExpToLevelSetxThreadStruct);

#endif /* _EXPTOLEVELSETX_H */
