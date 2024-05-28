/*!\file:  Chaco.h
 * \brief header file for Chaco module.
 */ 

#ifndef _CHACO_H
#define _CHACO_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
	#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*For python modules: needs to come before header files inclusion*/
#ifdef _HAVE_PYTHON_
#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#endif

/*headers*/
#include "../bindings.h" /*Should always come first to avoid python's warnings*/
#include <stdio.h>
#include <string.h>    /*  strcasecmp  */
#include <time.h>      /*  clock,time,difftime  */
#include "../../c/main/globals.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"

#undef __FUNCT__ 
#define __FUNCT__  "Chaco"

#ifdef _HAVE_MATLAB_MODULES_
/* serial input macros: */
#define A_IN      prhs[0]
#define VWGTS_IN  prhs[1]
#define EWGTS_IN  prhs[2]
#define X_IN      prhs[3]
#define Y_IN      prhs[4]
#define Z_IN      prhs[5]
#define OPTNS_IN  prhs[6]
#define NPARTS_IN prhs[7]
#define GOAL_IN   prhs[8]
/* serial output macros: */
#define ASSGN_OUT (mxArray**)&plhs[0]
#endif

#ifdef _HAVE_PYTHON_MODULES_
/* serial input macros: */
#define A_IN      PyTuple_GetItem(args,0)
#define VWGTS_IN  PyTuple_GetItem(args,1)
#define EWGTS_IN  PyTuple_GetItem(args,2)
#define X_IN      PyTuple_GetItem(args,3)
#define Y_IN      PyTuple_GetItem(args,4)
#define Z_IN      PyTuple_GetItem(args,5)
#define OPTNS_IN  PyTuple_GetItem(args,6)
#define NPARTS_IN PyTuple_GetItem(args,7)
#define GOAL_IN   PyTuple_GetItem(args,8)
/* serial output macros: */
#define ASSGN_OUT output,0
#endif

/* serial arg counts: */
#undef NLHS
#define NLHS  1
#undef NRHS
#define NRHS  9

#endif  /* _CHACO_H */
