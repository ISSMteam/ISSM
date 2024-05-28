/*
 * adolc_edf.h
 *
 *  Created on: Jun 26, 2012
 *      Author: utke
 */

#ifndef _ADOLC_EDF_H_
#define _ADOLC_EDF_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
#include "adolc/adolc.h"

struct Adolc_edf {
    ext_diff_fct *myEDF_for_solverx_p;
    Adolc_edf() : myEDF_for_solverx_p(0) {}
    inline friend std::ostream& operator << ( ostream&, const Adolc_edf& );
};

std::ostream& operator << ( std::ostream& out, const Adolc_edf& a) {
    out << a.myEDF_for_solverx_p;
    return out;
}

#endif

#endif
