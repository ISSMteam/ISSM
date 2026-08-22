/*! \file IntParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _INTPARAM_H_
#define _INTPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class IntParam: public Param{

	private:
		IssmInt value;

	public:
		/*IntParam constructors, destructors:*/
		IntParam();
		IntParam(int enum_type,IssmInt value);
		~IntParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return IntParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(int* pinteger){*pinteger=value;}
		void  SetValue(int integer){this->value=integer;}
};
#endif  /* _INTPARAM_H */
