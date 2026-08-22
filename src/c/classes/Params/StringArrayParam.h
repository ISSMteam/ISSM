/*! \file StringArrayParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _STRINGARRAYPARAM_H_
#define _STRINGARRAYPARAM_H_

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

class StringArrayParam: public Param{

	private: 
		char**   value;
		int      numstrings;

	public:
		/*StringArrayParam constructors, destructors: {{{*/
		StringArrayParam();
		StringArrayParam(int enum_type,char** values, int numstrings);
		~StringArrayParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(char*** pstringarray,int* pM);
		void  SetValue(char** stringarray,int M);
		/*}}}*/
};
#endif  /* _STRINGARRAYPARAM_H */
