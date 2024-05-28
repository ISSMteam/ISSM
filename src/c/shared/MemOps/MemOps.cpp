/*
 * MemOps.cpp
 *
 *  Created on: Sep 10, 2013
 *      Author: utke
 */

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "MemOps.h"

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> adouble*  xNew(unsigned int size, const char* const contig) {
	if (contig[0] == 't' || contig[0] == 'c')
		ensureContiguousLocations(size);

	adouble* aT_p=new adouble[size];
	assert(aT_p);
	return aT_p;
}
#endif
