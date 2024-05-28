#ifndef _BINDINGS_H_
#define _BINDINGS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef  _HAVE_PYTHON_MODULES_
#include "./python/include/pythonincludes.h"
#include "./python/include/wrapper_macros.h"
#include "./python/io/pythonio.h"
#endif

#ifdef  _HAVE_MATLAB_MODULES_
#include "./matlab/include/matlabincludes.h"
#include "./matlab/include/wrapper_macros.h"
#include "./matlab/io/matlabio.h"
#endif

#ifdef  _HAVE_JAVASCRIPT_MODULES_
#include "./javascript/include/javascriptincludes.h"
#include "./javascript/include/wrapper_macros.h"
#include "./javascript/io/javascriptio.h"
#endif


#endif
