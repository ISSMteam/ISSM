/*!\file InterpFromMesh2d.h
 * \brief: prototype for Data Interpolation mex module.
 */

#ifndef _INTERPFROMMESH2D_H
#define _INTERPFROMMESH2D_H

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
#define __FUNCT__  "InterpFromMesh2d"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEXHANDLE   prhs[0]
#define XHANDLE       prhs[1]
#define YHANDLE       prhs[2]
#define DATAHANDLE    prhs[3]
#define XPRIMEHANDLE  prhs[4]
#define YPRIMEHANDLE  prhs[5]
#define DEFAULTHANDLE prhs[6]
#define FILENAME      prhs[7]
/* serial output macros: */
#define DATAPRIME (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEXHANDLE   PyTuple_GetItem(args,0)
#define XHANDLE       PyTuple_GetItem(args,1)
#define YHANDLE       PyTuple_GetItem(args,2)
#define DATAHANDLE    PyTuple_GetItem(args,3)
#define XPRIMEHANDLE  PyTuple_GetItem(args,4)
#define YPRIMEHANDLE  PyTuple_GetItem(args,5)
#define DEFAULTHANDLE PyTuple_GetItem(args,6)
#define FILENAME      PyTuple_GetItem(args,7)
/* serial output macros: */
#define DATAPRIME output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1

#endif  /* _INTERPFROMMESH2D_H */
