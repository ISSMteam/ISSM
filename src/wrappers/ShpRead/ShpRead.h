/*!\file ShpRead.h
 * \brief: prototype for shp read mex module.
 */

#ifndef _SHPREAD_H
#define _SHPREAD_H

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
#define __FUNCT__  "ShpRead"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define SHP_IN  prhs[0]
/* serial output macros: */
#define SHP_OUT  (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define SHP_IN PyTuple_GetItem(args,0)
/* serial output macros: */
#define SHP_OUT output,0
#endif

/* serial arg counts: */
#undef NRHS
#define NRHS  1
#undef NLHS
#define NLHS  1

#endif
