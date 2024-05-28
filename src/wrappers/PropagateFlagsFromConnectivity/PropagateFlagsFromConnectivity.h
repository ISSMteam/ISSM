/*
	PropagateFlagsFromConnectivity.h
*/

#ifndef _PROPAGATEFLAGSFROMCONNECTIVITY_H
#define _PROPAGATEFLAGSFROMCONNECTIVITY_H

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
#define __FUNCT__  "PropagateFlagsFromConnectivity"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define CONNECTIVITY prhs[0]
#define POOL         prhs[1]
#define INDEX        prhs[2]
#define FLAGS        prhs[3]
/* serial output macros: */
#define POOLOUT (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define CONNECTIVITY PyTuple_GetItem(args,0)
#define POOL         PyTuple_GetItem(args,1)
#define INDEX        PyTuple_GetItem(args,2)
#define FLAGS        PyTuple_GetItem(args,3)
/* serial output macros: */
#define POOLOUT output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  4

#endif  /* _PROPAGATEFLAGSFROMCONNECTIVITY_H */
