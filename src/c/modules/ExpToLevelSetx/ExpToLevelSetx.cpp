/*! \file  ExpToLevelSetx.c
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ExpToLevelSetx.h"

int ExpToLevelSetx(double** pdistance,double* x, double* y, int nods, Contours* contours){

	/*output: */
	double* distance = xNew<double>(nods);
	for(int i=0;i<nods;i++) distance[i]=1e50;

	/*initialize thread parameters: */
	ExpToLevelSetxThreadStruct gate;
	gate.contours  = contours;
	gate.nods      = nods;
	gate.distance  = distance;
	gate.x         = x;
	gate.y         = y;

	/*launch the thread manager with ExpToLevelSetxt as a core: */
	LaunchThread(ExpToLevelSetxt,(void*)&gate,_NUMTHREADS_);

	/*Assign output pointers: */
	*pdistance=distance;

	return 1;
}
