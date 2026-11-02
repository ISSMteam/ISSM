/*! \file GPUHOParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _GPUHOPARAM_H_
#define _GPUHOPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Param.h"
#include "../../shared/shared.h"
class GPUHOParam: public Param{

	public:
		int Kff_dofarraysize;
		int Pf_valuesarraysize;
		int Kff_valuesarraysize;
		int PDStressHO_valuesarraysize;
		int PDStressHO_basisarraysize;
		int PDStressHO_factorarraysize;
		std::vector<int> fset_dofs;

		/*GPUHOParam constructors, destructors:*/
		GPUHOParam();
		~GPUHOParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return GPUHOParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(GPUHOParam** p_metadata){*p_metadata=this;};
};
#endif
