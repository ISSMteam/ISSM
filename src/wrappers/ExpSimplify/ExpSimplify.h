/*!\file ExpSimplify.h
 * \brief: prototype for exp to kml file conversion mex module.
 */

#ifndef _EXPSIMPLIFY_H
#define _EXPSIMPLIFY_H

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
#define __FUNCT__  "ExpSimplify"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define EXPFILE   prhs[0]
#define TOLERANCE prhs[1]
/* serial output macros: */
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define EXPFILE   PyTuple_GetItem(args,0)
#define TOLERANCE PyTuple_GetItem(args,1)
#endif

/* serial arg counts: */
#undef NRHS
#define NRHS  2
#undef NLHS
#define NLHS  0

#endif
