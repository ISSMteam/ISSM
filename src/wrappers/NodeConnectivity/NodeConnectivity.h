/*
	NodeConnectivity.h
*/

#ifndef _NODECONNECTIVITY_H
#define _NODECONNECTIVITY_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

#undef __FUNCT__ 
#define __FUNCT__  "NodeConnectivity"

/*Header files: */
#include "../bindings.h"
#include "../../c/main/globals.h"
#include "../../c/toolkits/toolkits.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../../c/shared/io/io.h"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define ELEMENTS prhs[0]
#define NUMNODES prhs[1]
/* serial output macros: */
#define CONNECTIVITY (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define ELEMENTS PyTuple_GetItem(args,0)
#define NUMNODES PyTuple_GetItem(args,1)
/* serial output macros: */
#define CONNECTIVITY output,0
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define ELEMENTS elementsin, nelsin,3
#define NUMNODES nodsin
/* serial output macros: */
#define CONNECTIVITY pconnectivity,pnods,pwidth
#define WRAPPER(modulename) extern "C" { int  NodeConnectivityModule(double** pconnectivity, int* pnods, int *pwidth, int* elementsin, int nelsin, int nodsin)
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  2

#endif  /* _NODECONNECTIVITY_H */
