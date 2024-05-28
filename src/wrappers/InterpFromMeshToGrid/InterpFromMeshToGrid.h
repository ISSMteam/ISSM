/*
	InterpFromMeshToGrid.h
*/

#ifndef _INTERPFROMMESHTOGRID_H
#define _INTERPFROMMESHTOGRID_H

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
#include "../../c/toolkits/toolkits.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../../c/shared/io/io.h"

#undef __FUNCT__ 
#define __FUNCT__  "InterpFromMeshToGrid"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEX        prhs[0]
#define X            prhs[1]
#define Y            prhs[2]
#define MESHDATA     prhs[3]
#define XGRID        prhs[4]
#define YGRID        prhs[5]
#define DEFAULTVALUE prhs[6]
/* serial output macros: */
#define GRIDDATA (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEX        PyTuple_GetItem(args,0)
#define X            PyTuple_GetItem(args,1)
#define Y            PyTuple_GetItem(args,2)
#define MESHDATA     PyTuple_GetItem(args,3)
#define XGRID        PyTuple_GetItem(args,4)
#define YGRID        PyTuple_GetItem(args,5)
#define DEFAULTVALUE PyTuple_GetItem(args,10)
/* serial output macros: */
#define GRIDDATA output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  7

#endif  /* _INTERPFROMMESHTOGRID_H*/
