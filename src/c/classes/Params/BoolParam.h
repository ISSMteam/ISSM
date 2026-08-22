/*! \file BoolParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _BOOLPARAM_H_
#define _BOOLPARAM_H_

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

class BoolParam: public Param{

	public:
		/*just hold 3 values for 3 vertices: */
		bool value;

		/*BoolParam constructors, destructors: {{{*/
		BoolParam();
		BoolParam(int enum_type,bool value);
		~BoolParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(bool* pbool){*pbool=value;}
		void  SetValue(bool boolean){this->value=boolean;}
		/*}}}*/
};
#endif  /* _BOOLPARAM_H */
