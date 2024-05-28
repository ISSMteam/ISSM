/*! \file DoubleTransientMatParam.h 
 *  \brief: header file for DoubleTransientMatParam object
 */

#ifndef _DOUBLETRANSIENTMATPARAM_H_
#define _DOUBLETRANSIENTMATPARAM_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"
/*}}}*/

class DoubleTransientMatParam: public DoubleMatParam{

	public:
		/*DoubleTransientMatParam constructors, destructors: {{{*/
		DoubleTransientMatParam(int enum_type,IssmDouble* value,int M,int N);
		/*}}}*/
};
#endif  /* _DOUBLETRANSIENTMATPARAM_H */
