/*! \file DataSetParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _DATASETPARAM_H_
#define _DATASETPARAM_H_

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

class DataSetParam: public Param{

	public:
		DataSet* value;

		/*DataSetParam constructors, destructors: {{{*/
		DataSetParam();
		DataSetParam(int enum_type,DataSet* dataset);
		~DataSetParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(DataSet** pdataset);
		void  SetValue(DataSet* dataset);
		/*}}}*/
};
#endif  /* _INTPARAM_H */
