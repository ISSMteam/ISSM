/*! \file DoubleVecParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _DOUBLEVECPARAM_H_
#define _DOUBLEVECPARAM_H_

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

class DoubleVecParam: public Param{

	private: 
		IssmDouble *values;
		int         M;

	public:
		/*DoubleVecParam constructors, destructors: {{{*/
		DoubleVecParam();
		DoubleVecParam(int enum_type,IssmDouble* values,int M);
		~DoubleVecParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum(){return DoubleVecParamEnum;}
		/*}}}*/
		/*Param virtual functions definitions: {{{*/
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, int* pN);
		void  SetValue(IssmDouble* IssmDoublearray,int M);
		/*}}}*/
		/*DoubleVecParam specific routines:{{{*/
		void  GetParameterValueByPointer(IssmDouble** pIssmDoublearray,int* pM);
		/*}}}*/

};
#endif  /* _DOUBLEVECPARAM_H */
