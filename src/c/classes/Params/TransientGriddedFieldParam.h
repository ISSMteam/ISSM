/*! \file TransientGriddedFieldParam.h 
 *  \brief: header file for transient gridded field 
 */

#ifndef _TRANSIENTGRIDDEDFIELDPARAM_H_
#define _TRANSIENTGRIDDEDFIELDPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class TransientGriddedFieldParam: public Param{

	protected:
		int         N;
		int         M;
		int         MN;
		int         T;
		bool        interpolation;
		bool        cycle;
		IssmDouble *values;
		IssmDouble *timesteps;

	public:
		/*TransientGriddedFieldParam: constructors, destructors:*/
		TransientGriddedFieldParam();
		TransientGriddedFieldParam(int in_enum_type,IssmDouble* in_values,IssmDouble* in_time,bool interpolation_on,bool cycle_in,int in_N,int in_M, int in_T);
		~TransientGriddedFieldParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return TransientGriddedFieldParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble time);
		void  GetParameterValue(IssmDouble* pdouble,int* index,int row,int column,IssmDouble time);
		void  GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble timestart,IssmDouble timeend);
		void  GetParameterValue(IssmDouble* pdouble,int row,int column,IssmDouble time, int timestepping, IssmDouble dt);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble time);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble starttime,IssmDouble endtime);
};
#endif  /* _TRANSIENTGRIDDEDFIELDPARAM_H_ */
