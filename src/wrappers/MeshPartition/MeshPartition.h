/*
	MeshPartition.h
*/

#ifndef _MESHPARTITION_H
#define _MESHPARTITION_H

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
#define __FUNCT__  "MeshPartition"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define NUMBEROFVERTICES   prhs[0]
#define ELEMENTS           prhs[1]
#define NUMBEROFVERTICES2D prhs[2]
#define ELEMENTS2D         prhs[3]
#define NUMBEROFLAYERS     prhs[4]
#define MESHELEMENTTYPE    prhs[5]
#define NUMAREAS           prhs[6]
/* serial output macros: */
#define ELEMENTPARTITIONING (mxArray**)&plhs[0]
#define NODEPARTITIONING (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define NUMBEROFVERTICES   PyTuple_GetItem(args,0)
#define ELEMENTS           PyTuple_GetItem(args,1)
#define NUMBEROFVERTICES2D PyTuple_GetItem(args,2)
#define ELEMENTS2D         PyTuple_GetItem(args,3)
#define NUMBEROFLAYERS     PyTuple_GetItem(args,4)
#define MESHELEMENTTYPE    PyTuple_GetItem(args,5)
#define NUMAREAS           PyTuple_GetItem(args,6)
/* serial output macros: */
#define ELEMENTPARTITIONING output,0
#define NODEPARTITIONING output,1
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS  7

#endif  /* _MESHPARTITION_H */
