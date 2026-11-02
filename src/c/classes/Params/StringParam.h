/*! \file StringParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _STRINGPARAM_H_
#define _STRINGPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class StringParam: public Param{

	private:
		/*just hold 3 values for 3 vertices: */
		char *value;

	public:
		/*StringParam constructors, destructors:*/
		StringParam();
		StringParam(int enum_type,char* value);
		~StringParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return StringParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(char** pstring);
		void  SetValue(char* string);
};
#endif  /* _STRINGPARAM_H */
