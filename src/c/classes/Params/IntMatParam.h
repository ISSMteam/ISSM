/*! \file IntMatParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _INTMATPARAM_H_
#define _INTMATPARAM_H_

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

class IntMatParam: public Param{

	private: 
		int* value;
		int M;
		int N;

	public:
		/*IntMatParam constructors, destructors: {{{*/
		IntMatParam();
		IntMatParam(int enum_type,int* value,int M,int N);
		~IntMatParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param vritual function definitions: {{{*/
		void  GetParameterValue(int** pintarray,int* pM,int* pN);
		void  SetValue(int* intarray,int M,int N);
		/*}}}*/
};
#endif  /* _INTMATPARAM_H */
