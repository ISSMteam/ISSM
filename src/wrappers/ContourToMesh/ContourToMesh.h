/*
	ContourToMesh.h
*/

#ifndef _CONTOURTOMESH_H
#define _CONTOURTOMESH_H

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
#define __FUNCT__ "ContourToMesh"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEX       prhs[0]
#define X           prhs[1]
#define Y           prhs[2]
#define CONTOUR     prhs[3]
#define INTERPTYPE  prhs[4]
#define EDGEVALUE   prhs[5]
/* serial output macros: */
#define PLHS0 (mxArray**)&plhs[0]
#define PLHS1 (mxArray**)&plhs[1]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEX       PyTuple_GetItem(args,0)
#define X           PyTuple_GetItem(args,1)
#define Y           PyTuple_GetItem(args,2)
#define CONTOUR     PyTuple_GetItem(args,3)
#define INTERPTYPE  PyTuple_GetItem(args,4)
#define EDGEVALUE   PyTuple_GetItem(args,5)
/* serial output macros: */
#define PLHS0 output,0
#define PLHS1 output,1
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define INDEX       indexin,nelin,3
#define X           xin,nodsin,1
#define Y           yin,nodsin,1
#define CONTOUR     contourx,contoury,contour_nods
#define INTERPTYPE  interptypein
#define EDGEVALUE   valuein
#define WRAPPER(modulename) extern "C" { int  ContourToMeshModule(double** pin_nod, double** pin_nel, double* indexin, double* xin, double* yin, double* contourx, double* contoury, char* interptypein, int nelin, int nodsin, int contour_nods, double valuein)
/* serial output macros: */
#define PLHS0 pin_nod,NULL
#define PLHS1 pin_nel,NULL
#define nrhs 6
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  2
#undef NRHS
#define NRHS 6

#endif  /* _CONTOURTOMESH_H */
