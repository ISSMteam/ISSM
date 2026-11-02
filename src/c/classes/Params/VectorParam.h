/*! \file VectorParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _VECTORPARAM_H_
#define _VECTORPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class VectorParam: public Param{

	private:
		Vector<IssmDouble>* value;

	public:
		/*VectorParam constructors, destructors:*/
		VectorParam();
		VectorParam(int enum_type,Vector<IssmDouble>* value);
		~VectorParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!"); };
		int    ObjectEnum(){return VectorParamEnum;}

		/*Param vritual function definitions:*/
		void  GetParameterValue(Vector<IssmDouble>** poutput);
		void  SetValue(Vector<IssmDouble>* vec);
};
#endif  /* _VECTORPARAM_H */
