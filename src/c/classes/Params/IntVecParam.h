/*! \file IntVecParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _INTVECPARAM_H_
#define _INTVECPARAM_H_

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

class IntVecParam: public Param{

	private: 
		int  M;
		int* values;

	public:
		/*IntVecParam constructors, destructors: {{{*/
		IntVecParam();
		IntVecParam(int enum_type,int* values,int M);
		IntVecParam(int enum_type,IssmDouble* values,int M);
		~IntVecParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual functions definitions: {{{*/
		void  GetParameterValue(int** pintarray,int* pM);
		void  SetValue(int* intarray,int M);
		/*}}}*/
};
#endif
