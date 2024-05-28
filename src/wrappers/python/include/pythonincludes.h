
#ifndef _PYTHON_INCLUDES_H_
#define _PYTHON_INCLUDES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_PYTHON_

#if _PYTHON_MAJOR_ >= 2
#undef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#else
#define NPY_NO_DEPRECATED_API 
#endif

#include <Python.h>
#include <arrayobject.h>

#endif
#endif /*_PYTHON_INCLUDES_H_*/
