/*!\file:  IssmConfig.h
 * \brief header file for IssmConfig module.
 */ 

#ifndef _ISSMCONFIG_H
#define _ISSMCONFIG_H

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
#include "../../c/shared/shared.h"

#undef __FUNCT__ 
#define __FUNCT__  "IssmConfig"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define NAME (mxArray*)prhs[0]
/* serial output macros: */
#define VALUE (mxArray**)&plhs[0]
#define SVALUE (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define NAME PyTuple_GetItem(args,0)
/* serial output macros: */
#define VALUE output,0
#define SVALUE output,0
#endif


#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define NAME string
/* serial output macros: */
#define VALUE pvalue
#define SVALUE psvalue
#define WRAPPER(modulename) extern "C" { int IssmConfigModule(double* pvalue, char** psvalue, char* string)
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  1

#endif  /* _TEST_H */
