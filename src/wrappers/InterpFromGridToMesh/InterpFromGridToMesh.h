/*!\file InterpFromGridToMesh.h
 * \brief: prototype for Data Interpolation mex module.
 */

#ifndef _InterpFromGridToMesh_H
#define _InterpFromGridToMesh_H

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
#define __FUNCT__  "InterpFromGridToMesh"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define XHANDLE       prhs[0]
#define YHANDLE       prhs[1]
#define DATAHANDLE    prhs[2]
#define XMESHHANDLE   prhs[3]
#define YMESHHANDLE   prhs[4]
#define DEFAULTHANDLE prhs[5]
#define INTERPENUM    prhs[6]
/* serial output macros: */
#define DATAMESH (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define XHANDLE       PyTuple_GetItem(args,0)
#define YHANDLE       PyTuple_GetItem(args,1)
#define DATAHANDLE    PyTuple_GetItem(args,2)
#define XMESHHANDLE   PyTuple_GetItem(args,3)
#define YMESHHANDLE   PyTuple_GetItem(args,4)
#define DEFAULTHANDLE PyTuple_GetItem(args,5)
#define INTERPENUM    PyTuple_GetItem(args,6)
/* serial output macros: */
#define DATAMESH output,0
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define XHANDLE       xIn,dataNumColsIn,1
#define YHANDLE       yIn,dataNumRowsIn,1
#define DATAHANDLE    dataIn,dataNumRowsIn,dataNumColsIn
#define XMESHHANDLE   xMeshIn,meshNumRowsIn,1
#define YMESHHANDLE   yMeshIn,meshNumRowsIn,1
#define DEFAULTHANDLE defaultValue
#define INTERPENUM    interpType
/* serial output macros: */
#define DATAMESH pdataMesh
#define WRAPPER(modulename) extern "C" { int InterpFromGridToMeshModule(double** pdataMesh, double* xIn, double* yIn, double* dataIn, double* xMeshIn, double* yMeshIn, double defaultValue, int nodsIn, int dataNumRowsIn, int dataNumColsIn, int meshNumRowsIn, char* interpType)
#define nrhs 6
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  6

#endif  /* _INTERPFROMGRIDTOMESH_H */
