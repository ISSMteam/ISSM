/*
	CoordTransform.h
*/

#ifndef _COORDTRANSFORM_H
#define _COORDTRANSFORM_H

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
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"

#undef __FUNCT__ 
#define __FUNCT__  "CoordTransform"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define XIN     prhs[0]
#define YIN     prhs[1]
#define PROJIN  prhs[2]
#define PROJOUT prhs[3]
/* serial output macros: */
#define XOUT (mxArray**)&plhs[0]
#define YOUT (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define XIN     PyTuple_GetItem(args,0)
#define YIN     PyTuple_GetItem(args,1)
#define PROJIN  PyTuple_GetItem(args,2)
#define PROJOUT PyTuple_GetItem(args,3)
/* serial output macros: */
#define XOUT output,0
#define YOUT output,1
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define XIN     xin
#define YIN     yin
#define PROJIN  projin
#define PROJOUT projout
/* serial output macros: */
/*NOT IMPLEMENTED YET*/
#endif


/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS  4

#endif
