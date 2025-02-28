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
		IssmDouble* gradient; 
		int         enum_type;
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
		void  GetParameterValue(bool* pbool){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a bool");}
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
		int   InstanceEnum(){return enum_type;}

		void  SetEnum(int enum_in){this->enum_type = enum_in;};
		void  SetValue(bool boolean){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a boolean");}
		void  SetValue(int integer){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an integer");}
		void  SetValue(IssmDouble scalar){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a scalar");}
		void  SetValue(char* string){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a string");}
		void  SetValue(char** stringarray,int M){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a string array");}
		void  SetValue(IssmDouble* IssmDoublearray){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a scalar");}
		void  SetValue(IssmDouble* IssmDoublearray,int M){_error_("not implemented");}
		void  SetValue(IssmDouble* pIssmDoublearray,int M,int N);
		void  SetValue(int* intarray,int M){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an integer");}
		void  SetValue(int* pintarray,int M,int N){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a int mat array");}
		void  SetValue(Vector<IssmDouble>* vec){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a Vec");}
		void  SetValue(Matrix<IssmDouble>* mat){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a Mat");}
		void  SetValue(FILE* fid){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a FILE");}
		void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an array of matrices");}
		void  SetGradient(IssmDouble* poutput, int M, int N);
		void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int N,const char* data,int offset);
		/*}}}*/
};
#endif  /* _DOUBLEPARAM_H */
