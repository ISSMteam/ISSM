/*! \file DoubleMatParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _DOUBLEMATPARAM_H_
#define _DOUBLEMATPARAM_H_

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

class DoubleMatParam: public Param{

	protected: 
		IssmDouble *value;
		int         M;
		int         N;

	public:
		/*DoubleMatParam constructors, destructors: {{{*/
		DoubleMatParam();
		DoubleMatParam(int enum_type,IssmDouble* value,int M,int N);
		~DoubleMatParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		void  Echo();
		void  DeepEcho();
		int   ObjectEnum();
		Param* copy();
		void Marshall(MarshallHandle* marshallhandle);
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN);

		void  SetValue(IssmDouble* IssmDoublearray,int M,int N);
		/*}}}*/
		/*DoubleMatParam specific routines:{{{*/
		void  GetParameterValueByPointer(IssmDouble** pIssmDoublearray,int* pM,int* pN);
		/*}}}*/
};
#endif  /* _DOUBLEMATPARAM_H */
