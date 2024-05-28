/*
	ElementConnectivity.h
*/

#ifndef _ELEMENTCONNECTIVITY_H
#define _ELEMENTCONNECTIVITY_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

/*Header files: */
#include "../bindings.h"
#include "../../c/main/globals.h"
#include "../../c/toolkits/toolkits.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../../c/shared/io/io.h"

#undef __FUNCT__ 
#define __FUNCT__  "ElementConnectivity"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define ELEMENTS         prhs[0]
#define NODECONNECTIVITY prhs[1]
/* serial output macros: */
#define ELEMENTCONNECTIVITY (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define ELEMENTS         PyTuple_GetItem(args,0)
#define NODECONNECTIVITY PyTuple_GetItem(args,1)
/* serial output macros: */
#define ELEMENTCONNECTIVITY output,0
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define ELEMENTS         elementsin, nelsin,3
#define NODECONNECTIVITY nodeconnectivityin, nodsin, widthin
/* serial output macros: */
#define ELEMENTCONNECTIVITY pelementconnectivity,NULL,NULL
#define WRAPPER(modulename) extern "C" { int  ElementConnectivityModule(double** pelementconnectivity, int* elementsin, int* nodeconnectivityin, int nelsin, int nodsin, int widthin)
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  2

#endif  /* _ELEMENTCONNECTIVITY_H */
