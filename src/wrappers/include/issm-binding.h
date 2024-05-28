#ifndef _ISSM_BINDING_H_
#define _ISSM_BINDING_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef  _HAVE_MATLAB_MODULES_
#include "../matlab/include/matlab-macros.h"
#endif

#ifdef  _HAVE_PYTHON_MODULES_
#include "../python/include/python-macros.h"
#endif

#endif
