/*! \file ControlParam.h 
 *  \brief: header file for ControlParam object
 */

#ifndef _CONTROLPARAM_H_
#define _CONTROLPARAM_H_

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

class ControlParam: public Param{

	private: 
		IssmDouble* value;   //Can either be a VecParam or a TransientParam
		IssmDouble* minvalue;
		IssmDouble* maxvalue;
		int         M,N;

	public:
		/*ControlParam constructors, destructors: {{{*/
		ControlParam();
		ControlParam(IssmDouble* in_value, IssmDouble* in_minvalue, IssmDouble* in_maxvalue, int in_enum_type, int in_M, int in_N);
		~ControlParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual functions definitions: {{{*/
		void  GetParameterValue(int* pinteger){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an integer");}
		void  GetParameterValue(int** pintarray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an integer");}
		void  GetParameterValue(int** pintarray,int* pM,int* pN){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a matrix");}
		void  GetParameterValue(IssmDouble* pIssmDouble);
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time);
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){_error_("Param "<< EnumToStringx(enum_type) << " does not support cycling of forcing");}
		void  GetParameterValue(char** pstring){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string");}
		void  GetParameterValue(char*** pstringarray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string array");}
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, int* pN){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble array");}
		void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a matrix array");}
		void  GetParameterValue(Vector<IssmDouble>** pvec){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a Vec");}
		void  GetParameterValue(Matrix<IssmDouble>** pmat){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a Mat");}
		void  GetParameterValue(FILE** pfid){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a FILE");}
		void  GetParameterValue(DataSet** pdataset){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a DataSet");}
		void  SetValue(IssmDouble* pIssmDoublearray,int M,int N);
		void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int N,const char* data,int offset);
		/*}}}*/
};
#endif  /* _DOUBLEPARAM_H */
