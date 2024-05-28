/* \file matlab macros.h
 * \brief: macros used for the matlab bindings
 */

#ifndef _MATLAB_MACROS_H_
#define _MATLAB_MACROS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_MATLAB_
/* MODULEBOOT/MODULEEND {{{*/

/*The following macros hide the error exception handling in a matlab module. Just put 
 * MODULEBOOT(); and MODULEEND(); at the beginning and end of a module, and c++ exceptions 
 * will be trapped*/
#define MODULEBOOT(); try{ \
	IssmComm::SetComm(); \
	ToolkitOptions::Init();

#define MODULEEND(); }\
	catch(ErrorException &exception){\
		mexErrMsgTxt(exception.WrapperReport()); \
	}\
	catch (exception &e){\
		mexErrMsgTxt(e.what());\
	}\
	catch(...){\
		mexErrMsgTxt("An unexpected error occurred");\
	}
/*}}}*/
/* WRAPPER {{{*/
#define WRAPPER(modulename,...) void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) 
/*}}}*/
/* CHECKARGUMENTS {{{*/
#define CHECKARGUMENTS(NLHS,NRHS,functionpointer) CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,functionpointer)
/*}}}*/
#endif

#endif
