/* \file python_macros.h
 * \brief: macros used for the python bindings
 */

#ifndef _PY_WRAPPER_MACROS_H_
#define _PY_WRAPPER_MACROS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_PYTHON_
/* MODULEBOOT/MODULEEND {{{*/

/*The following macros hide the error exception handling in a python module. Just put 
 * MODULEBOOT(); and MODULEEND(); at the beginning and end of a module, and c++ exceptions 
 * will be trapped*/
#define MODULEBOOT(); \
	PyObject *output = PyTuple_New(NLHS); \
	int       nrhs   = (int)PyTuple_Size(args);  \
	if(!output) return NULL;\
	try{ \
	IssmComm::SetComm();

#define MODULEEND(); }\
  catch(ErrorException &exception){\
	  PyErr_SetString(PyExc_TypeError,exception.WrapperReport()); \
	  return NULL;\
  } \
	catch (exception &e){\
		PyErr_SetString(PyExc_TypeError,e.what());\
		return NULL;\
	}\
	catch(...){\
		PyErr_SetString(PyExc_TypeError,"An unexpected error occurred");\
		return NULL;\
	}\
	return output;
//}}}
#if _PYTHON_MAJOR_ >=3
/* WRAPPER 3.2 {{{*/
#define WRAPPER(modulename,...)  \
static PyObject* modulename(PyObject* self,PyObject* args);\
static PyMethodDef modulename##_funcs[] = {\
	{#modulename, (PyCFunction)modulename, METH_VARARGS, ""},\
	{NULL,NULL,0,NULL}\
};\
static struct PyModuleDef modulename##module= {\
	PyModuleDef_HEAD_INIT,\
	#modulename,   /* name of module */\
	NULL, /* module documentation, may be NULL */\
	-1,       /* size of per-interpreter state of the module,\
				 or -1 if the module keeps state in global variables. */\
	modulename##_funcs\
};\
PyMODINIT_FUNC PyInit_##modulename(void){\
	import_array();\
	return PyModule_Create(&modulename##module);\
}\
static PyObject* modulename(PyObject* self,PyObject* args)
/*}}}*/
#else
/* WRAPPER 2.7 {{{*/
#define WRAPPER(modulename,...)  \
static PyObject* modulename(PyObject* self,PyObject* args);\
static PyMethodDef modulename##_funcs[] = {\
	{#modulename, (PyCFunction)modulename, METH_VARARGS, ""},\
	{NULL,NULL,0,NULL}\
};\
PyMODINIT_FUNC init##modulename(void){\
	import_array();\
	(void) Py_InitModule(#modulename, modulename##_funcs);\
}\
static PyObject* modulename(PyObject* self,PyObject* args)
/*}}}*/
#endif
/* CHECKARGUMENTS {{{*/
#define CHECKARGUMENTS(NLHS,NRHS,functionpointer) CheckNumPythonArguments(args, NRHS,functionpointer)
/*}}}*/
#endif
#endif
