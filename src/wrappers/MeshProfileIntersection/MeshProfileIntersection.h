/*
	MeshProfileIntersection.h
*/

#ifndef _MESHPROFILEINTERSECTION_H
#define _MESHPROFILEINTERSECTION_H

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
#define __FUNCT__ "MeshProfileIntersection"

#ifdef _HAVE_MATLAB_MODULES_
/* input macros: */
#define INDEX    prhs[0]
#define X        prhs[1]
#define Y        prhs[2]
#define FILENAME prhs[3]
/* serial output macros: */
#define SEGMENTS (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* input macros: */
#define INDEX    PyTuple_GetItem(args,0)
#define X        PyTuple_GetItem(args,1)
#define Y        PyTuple_GetItem(args,2)
#define FILENAME PyTuple_GetItem(args,3)
/* serial output macros: */
#define SEGMENTS output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS 1
#undef NRHS
#define NRHS 4

#endif  /* _MESHPROFILEINTERSECTION_H */
