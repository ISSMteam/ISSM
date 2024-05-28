/*
	ExpToLevelSet.h
*/

#ifndef _EXPTOLEVELSET_H
#define _EXPTOLEVELSET_H

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
#define __FUNCT__ "ExpToLevelSet"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define X           prhs[0]
#define Y           prhs[1]
#define CONTOUR     prhs[2]

/* serial output macros: */
#define PLHS0 (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define X           PyTuple_GetItem(args,0)
#define Y           PyTuple_GetItem(args,1)
#define CONTOUR     PyTuple_GetItem(args,2)
/* serial output macros: */
#define PLHS0 output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS 3

#endif  /* _EXPTOLEVELSET_H */
