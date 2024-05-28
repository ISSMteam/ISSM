/*!\file:  Scotch.h
 * \brief header file for Scotch module.
 */ 

#ifndef _SCOTCH_H
#define _SCOTCH_H

#include "../bindings.h" /*Should always come first to avoid python's warnings*/
#include <stdio.h>
#include <string.h>    /*  strcasecmp  */
#include <time.h>      /*  clock,time,difftime  */
#include "../../c/main/globals.h"
#include "../../c/modules/modules.h"
#include "../../c/shared/shared.h"

#undef __FUNCT__ 
#define __FUNCT__  "Scotch"

/*  Scotch structures and prototypes  */
#ifdef MATLAB
#include "mat.h"
#include "mex.h"
#include "matrix.h"

#define printf mexPrintf
#define fprintf(file,...) (file == stdout || file == stderr ? mexPrintf(__VA_ARGS__) : fprintf(file,__VA_ARGS__))
#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree
#define exit(status) mexErrMsgTxt("exit=" #status)
#endif

#endif  /* _SCOTCH_H */
