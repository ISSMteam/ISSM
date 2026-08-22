/*! \file MatrixParam.h 
 *  \brief: header file for MatrixParam object
 */

#ifndef _MATRIXPARAM_H_
#define _MATRIXPARAM_H_

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

class MatrixParam: public Param{

	private: 
		Matrix<IssmDouble>* value;

	public:
		/*MatrixParam constructors, destructors: {{{*/
		MatrixParam();
		MatrixParam(int enum_type,Matrix<IssmDouble>* value);
		~MatrixParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!"); };
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(Matrix<IssmDouble>** poutput);
		void  SetValue(Matrix<IssmDouble>* mat);
		/*}}}*/
};
#endif  /* _MATRIXPARAM_H */
