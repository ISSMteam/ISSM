/*
	Triangle.h
*/

#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

#if _IS_MAC_ == 0
#ifdef _HAVE_JAVASCRIPT_MODULES_
#undef _DO_NOT_LOAD_GLOBALS_ /*only module where this needs to be undefined, so as to 
							   not include IssmComm several times in the JavaScript module construct.*/
#endif
#endif

/*Header files: */
#include "../bindings.h"
#include "../../c/main/globals.h"
#include "../../c/toolkits/toolkits.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"
#include "../../c/shared/io/io.h"

#undef __FUNCT__ 
#define __FUNCT__  "Triangle"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define DOMAINOUTLINE  prhs[0]
#define RIFTSOUTLINE   prhs[1]
#define AREA           prhs[2]
/* serial output macros: */
#define INDEX             (mxArray**)&plhs[0]
#define X                 (mxArray**)&plhs[1]
#define Y                 (mxArray**)&plhs[2]
#define SEGMENTS          (mxArray**)&plhs[3]
#define SEGMENTMARKERLIST (mxArray**)&plhs[4]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define DOMAINOUTLINE PyTuple_GetItem(args,0)
#define RIFTSOUTLINE  PyTuple_GetItem(args,1)
#define AREA          PyTuple_GetItem(args,2)
/* serial output macros: */
#define INDEX             output,0
#define X                 output,1
#define Y                 output,2
#define SEGMENTS          output,3
#define SEGMENTMARKERLIST output,4
#endif

#ifdef _HAVE_JAVASCRIPT_MODULES_
/* serial input macros: */
#define DOMAINOUTLINE domainx,domainy,domainnods
#define RIFTSOUTLINE  NULL,NULL,0
#define AREA          areain
/* serial output macros: */
#define INDEX             pindex,pnel
#define X                 px,pnods
#define Y                 py,pnods
#define SEGMENTS          psegments,pnsegs
#define SEGMENTMARKERLIST psegmentmarkers,pnsegs
#define WRAPPER(modulename) extern "C" { int  TriangleModule(double** pindex, double** px, double** py, int* pnel, int* pnods, double** psegments, double** psegmentmarkers, int* pnsegs, double* domainx, double* domainy, int domainnods, double areain)
#define _DO_NOT_LOAD_GLOBALS_//we only load globals for TriangleModule.js, not other modules!
#endif


/* serial arg counts: */
#undef NLHS
#define NLHS  5
#undef NRHS
#define NRHS  3

#endif  /* _TRIANGLE_H */
