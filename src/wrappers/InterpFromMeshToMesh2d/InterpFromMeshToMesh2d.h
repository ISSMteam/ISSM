/*!\file InterpFromMeshToMesh2d.h
 * \brief: prototype for Data Interpolation mex module.
 */

#ifndef _INTERPFROMMESHTOMESH2d_H
#define _INTERPFROMMESHTOMESH2d_H

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
#define __FUNCT__  "InterpFromMeshToMesh2d"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEX     prhs[0]
#define X         prhs[1]
#define Y         prhs[2]
#define DATA      prhs[3]
#define XINTERP   prhs[4]
#define YINTERP   prhs[5]
#define ARGUMENTS prhs 
/* serial output macros: */
#define DATAINTERP (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEX          PyTuple_GetItem(args,0)
#define X              PyTuple_GetItem(args,1)
#define Y              PyTuple_GetItem(args,2)
#define DATA           PyTuple_GetItem(args,3)
#define XINTERP        PyTuple_GetItem(args,4)
#define YINTERP        PyTuple_GetItem(args,5)
#define ARGUMENTS args
/* serial output macros: */
#define DATAINTERP output,0
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define INDEX          indexin,nelin,3
#define X              xin,nodsin
#define Y              yin,nodsin
#define DATA           datain,nodsin,1
#define XINTERP        x_interpin, nods_interpin
#define YINTERP        y_interpin, nods_interpin
#define ARGUMENTS "default_value",default_value
/* serial output macros: */
#define DATAINTERP pdata_interp,NULL,NULL
#define WRAPPER(modulename) extern "C" { int  InterpFromMeshToMesh2dModule(double** pdata_interp,int* indexin,double* xin,double* yin,double* datain,double* x_interpin,double* y_interpin,int nelin,int nodsin,int nods_interpin,double default_value)
#define nrhs  6
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  6

#endif
