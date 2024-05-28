/*
	Kriging.h
*/

#ifndef _KRIGING_H_
#define _KRIGING_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

#include "../bindings.h"
#include "../../c/main/globals.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"

#undef __FUNCT__ 
#define __FUNCT__  "Kriging"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define X            prhs[0]
#define Y            prhs[1]
#define OBSERVATIONS prhs[2]
#define XINTERP      prhs[3]
#define YINTERP      prhs[4]

/* serial output macros: */
#define PREDICTIONS (mxArray**)&plhs[0]
#define ERROR       (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define X            PyTuple_GetItem(args,0)
#define Y            PyTuple_GetItem(args,1)
#define OBSERVATIONS PyTuple_GetItem(args,2)
#define XINTERP      PyTuple_GetItem(args,3)
#define YINTERP      PyTuple_GetItem(args,4)

/* serial output macros: */
#define PREDICTIONS output,0
#define ERROR       output,1
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS  5

#endif  /* _KRIGING_H_ */
