/*! \file  DistanceToMaskBoundaryx.c
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./DistanceToMaskBoundaryx.h"

int DistanceToMaskBoundaryx(IssmDouble** pdistance,IssmDouble* x, IssmDouble* y, IssmDouble* mask, int nods) {

	/*output: */
	IssmDouble*  distance;

	/*initialize: */
	distance=xNew<IssmDouble>(nods);

	/*initialize thread parameters: */
	DistanceToMaskBoundaryxThreadStruct gate;
	gate.distance = distance;
	gate.x        = x;
	gate.y        = y;
	gate.mask     = mask;
	gate.nods     = nods;

	/*launch the thread manager with DistanceToMaskBoundaryxt as a core: */
	LaunchThread(DistanceToMaskBoundaryxt,(void*)&gate,_NUMTHREADS_);

	/*Assign output pointers: */
	*pdistance=distance;

	return 1;
}
