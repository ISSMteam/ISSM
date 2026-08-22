/*! \file DoubleMatArrayParam.h 
 *  \brief: header file for object holding an array of serial matrices
 */

#ifndef _DOUBLEMATARRAYPARAM_H_
#define _DOUBLEMATARRAYPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class DoubleMatArrayParam: public Param{

	private:
		IssmDouble **array;        //array of matrices
		int          M;            //size of array
		int         *mdim_array;   //m-dimensions of matrices in the array
		int         *ndim_array;   //n-dimensions -f matrices in the array

	public:
		/*DoubleMatArrayParam constructors, destructors:*/
		DoubleMatArrayParam();
		DoubleMatArrayParam(int enum_type,IssmDouble** array, int M, int* mdim_array, int* ndim_array);
		~DoubleMatArrayParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return DoubleMatArrayParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims);
		void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array);
};
#endif  /* _DOUBLEMATARRAYPARAM_H */
