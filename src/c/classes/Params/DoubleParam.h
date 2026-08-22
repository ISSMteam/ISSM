/*! \file DoubleParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _DOUBLEPARAM_H_
#define _DOUBLEPARAM_H_

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

class DoubleParam: public Param{

	private:
		/*just hold 3 values for 3 vertices: */
		IssmDouble value;

	public:
		/*DoubleParam constructors, destructors: {{{*/
		DoubleParam();
		DoubleParam(int enum_type,IssmDouble value);
		~DoubleParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(IssmDouble* pIssmDouble){*pIssmDouble=value;};
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){*pdouble=value;};
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){*pdouble=value;};
		void  SetValue(bool boolean){this->value=(IssmDouble)boolean;}
		void  SetValue(int integer){this->value=(IssmDouble)integer;}
		void  SetValue(IssmDouble scalar){this->value=(IssmDouble)scalar;}
		/*}}}*/
};
#endif  /* _DOUBLEPARAM_H */
