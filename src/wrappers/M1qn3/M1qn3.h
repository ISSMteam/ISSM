/*!\file M1qn3.h
 * \brief: prototype for Data Interpolation mex module.
 */

#ifndef _M1QN3_WRAPPER_H
#define _M1QN3_WRAPPER_H

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
#define __FUNCT__  "M1qn3"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define XHANDLE       prhs[0]
#define GHANDLE       prhs[1]
#define JHANDLE       prhs[2]
/* serial output macros: */
#define XOUT (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define XHANDLE   PyTuple_GetItem(args,0)
#define GHANDLE   PyTuple_GetItem(args,1)
#define JHANDLE   PyTuple_GetItem(args,2)
/* serial output macros: */
#define XOUT output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1

#endif  /* _M1QN3_WRAPPER_H */
