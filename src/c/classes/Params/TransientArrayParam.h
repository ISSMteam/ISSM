/*! \file TransientArrayParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _TRANSIENTARRAYPARAM_H_
#define _TRANSIENTARRAYPARAM_H_

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

class TransientArrayParam: public Param{

	protected: 
		int         N;
		int         M;
		bool        interpolation;
		bool        cycle;
		IssmDouble *values;
		IssmDouble *timesteps;

	public:
		/*TransientArrayParam constructors, destructors: {{{*/
		TransientArrayParam();
		TransientArrayParam(int in_enum_type,IssmDouble* in_values,IssmDouble* in_time,bool interpolation_on,bool cycle_in,int in_N,int in_M);
		~TransientArrayParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum(){return TransientArrayParamEnum;}
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(IssmDouble* pdouble,int row,IssmDouble time);
		void  GetParameterValue(IssmDouble* pdouble,int row,IssmDouble time, int timestepping, IssmDouble dt);
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){_error_("Parameter " <<EnumToStringx(enum_type) << " needs row to be specified");}
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time,int timestepping, IssmDouble dt){_error_("Parameter " <<EnumToStringx(enum_type) << " needs row to be specified");}
		/*}}}*/
};
#endif  /* _TRANSIENTARRAYPARAM_H */
