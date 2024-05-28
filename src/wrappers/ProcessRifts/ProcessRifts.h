/*
 * ProcessRifts.h
 */ 

#ifndef _PROCESSRIFTS_H_
#define _PROCESSRIFTS_H_

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
#define __FUNCT__  "ProcessRifts"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define INDEXIN          prhs[0]
#define XIN              prhs[1]
#define YIN              prhs[2]
#define SEGMENTSIN       prhs[3]
#define SEGMENTMARKERSIN prhs[4]
/* serial output macros: */
#define INDEXOUT          (mxArray**)&plhs[0]
#define XOUT              (mxArray**)&plhs[1]
#define YOUT              (mxArray**)&plhs[2]
#define SEGMENTSOUT       (mxArray**)&plhs[3]
#define SEGMENTMARKERSOUT (mxArray**)&plhs[4]
#define RIFTSTRUCT        (mxArray**)&plhs[5]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define INDEXIN          PyTuple_GetItem(args,0)
#define XIN              PyTuple_GetItem(args,1)
#define YIN              PyTuple_GetItem(args,2)
#define SEGMENTSIN       PyTuple_GetItem(args,3)
#define SEGMENTMARKERSIN PyTuple_GetItem(args,4)
/* serial output macros: */
#define INDEXOUT          output,0
#define XOUT              output,1
#define YOUT              output,2
#define SEGMENTSOUT       output,3
#define SEGMENTMARKERSOUT output,4
#define RIFTSTRUCT        output,5
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  6
#undef NRHS
#define NRHS  5

#endif
