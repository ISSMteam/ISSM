/* \file javascript macros.h
 * \brief: macros used for the javascript bindings
 */

#ifndef _JAVASCRIPT_MACROS_H_
#define _JAVASCRIPT_MACROS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_JAVASCRIPT_
/* MODULEBOOT/MODULEEND {{{*/

/*The following macros hide the error exception handling in a javascript module. Just put 
 * MODULEBOOT(); and MODULEEND(); at the beginning and end of a module, and c++ exceptions 
 * will be trapped*/
#define MODULEBOOT(); try{ \
	IssmComm::SetComm();

#define MODULEEND(); }\
	catch(ErrorException &exception){\
		printf(exception.WrapperReport()); \
	}\
	catch (exception &e){\
		printf(e.what());\
	}\
	catch(...){\
		printf("An unexpected error occurred");\
	}\
	return 0;\
	}
/*}}}*/
/* CHECKARGUMENTS {{{*/
#define CHECKARGUMENTS(NLHS,NRHS,functionpointer)  //do nothing, we are not creating a dynamic library here!
/*}}}*/
#endif

#endif
