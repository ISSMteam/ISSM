/*! \file TransientParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _TRANSIENTPARAM_H_
#define _TRANSIENTPARAM_H_

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

class TransientParam: public Param{

	protected: 
		int         N;
		bool        interpolation;
		bool        cycle;
		IssmDouble *values;
		IssmDouble *timesteps;

	public:
		/*TransientParam constructors, destructors: {{{*/
		TransientParam();
		TransientParam(int in_enum_type,IssmDouble* in_values,IssmDouble* in_time,bool interpolation_in,bool cycle_in,int in_N);
		~TransientParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time);
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time,int timestepping,IssmDouble dt);
		/*}}}*/
};
#endif  /* _TRANSIENTPARAM_H */
