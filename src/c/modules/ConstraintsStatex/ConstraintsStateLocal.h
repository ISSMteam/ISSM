/*!\file:  ConstraintsStateLocal.h
 * \brief local header files
 */ 

#ifndef _CONSTRAINTSSTATELOCAL_H
#define _CONSTRAINTSSTATELOCAL_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../classes/classes.h"

/*rifts module: */
void RiftConstraintsState(int* pconverged, int* pnum_unstable_constraints,Loads* loads,int min_mechanical_constraints,int analysis_type);
void RiftConstrain(int* pnum_unstable_constraints,Loads* loads,int analysis_type);
int  RiftIsFrozen(Loads* loads,int analysis_type);
void RiftFreezeConstraints(Loads* loads,int analysis_type);

#endif  /* _CONSTRAINTSSTATEX_H */
