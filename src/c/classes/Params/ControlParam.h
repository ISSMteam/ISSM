/*! \file ControlParam.h 
 *  \brief: header file for ControlParam object
 */

#ifndef _CONTROLPARAM_H_
#define _CONTROLPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class ControlParam: public Param{

	private:
		IssmDouble* value;
		IssmDouble* minvalue;
		IssmDouble* maxvalue;
		int         M,N;

	public:
		/*ControlParam constructors, destructors:*/
		ControlParam();
		ControlParam(IssmDouble* in_value, IssmDouble* in_minvalue, IssmDouble* in_maxvalue, int in_enum_type, int in_M, int in_N);
		~ControlParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return ControlParamEnum;}

		/*Param virtual functions definitions:*/
		void  GetParameterValue(IssmDouble* pIssmDouble);
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data);
		void  SetValue(IssmDouble* pIssmDoublearray,int M,int N);
		void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int N,const char* data,int offset);
};
#endif  /* _DOUBLEPARAM_H */
