/*!\file BamgConvertMesh.h
 * \brief: prototype for Data Interpolation mex module.
 */

#ifndef _BAMGCONVERTMESH_H
#define _BAMGCONVERTMESH_H

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
#include "../../c/shared/io/io.h"

#undef __FUNCT__ 
#define __FUNCT__  "BamgConvertMesh"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEXHANDLE prhs[0]
#define XHANDLE     prhs[1]
#define YHANDLE     prhs[2]
/* serial output macros: */
#define BAMGMESHOUT    (mxArray**)&plhs[0]
#define BAMGGEOMOUT    (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEXHANDLE PyTuple_GetItem(args,0)
#define XHANDLE     PyTuple_GetItem(args,1)
#define YHANDLE     PyTuple_GetItem(args,2)
/* serial output macros: */
#define BAMGMESHOUT    output,0
#define BAMGGEOMOUT    output,1
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS  3

#endif
