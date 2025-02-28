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
		int        enum_type;
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
		int   Id(); 
		void Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		void  GetParameterValue(bool* pbool);
		void  GetParameterValue(int* pinteger);
		void  GetParameterValue(int** pintarray,int* pM);
		void  GetParameterValue(int** pintarray,int* pM,int* pN);
		void  GetParameterValue(IssmDouble* pIssmDouble){*pIssmDouble=value;};
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){*pdouble=value;};
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){*pdouble=value;};
		void  GetParameterValue(char** pstring){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string");}
		void  GetParameterValue(char*** pstringarray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string array");}
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, int* pN);
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble array");}
		void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a matrix array");}
		void  GetParameterValue(Vector<IssmDouble>** pvec){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a Vec");}
		void  GetParameterValue(Matrix<IssmDouble>** pmat){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a Mat");}
		void  GetParameterValue(FILE** pfid){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a FILE");}
		void  GetParameterValue(DataSet** pdataset){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a DataSet");}
		int   InstanceEnum(){return enum_type;}

		void  SetEnum(int enum_in){this->enum_type = enum_in;};
		void  SetValue(bool boolean){this->value=(IssmDouble)boolean;}
		void  SetValue(int integer){this->value=(IssmDouble)integer;}
		void  SetValue(IssmDouble scalar){this->value=(IssmDouble)scalar;}
		void  SetValue(char* string){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a string");}
		void  SetValue(char** stringarray,int M){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a string array");}
		void  SetValue(IssmDouble* IssmDoublearray){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a IssmDouble array");}
		void  SetValue(IssmDouble* IssmDoublearray,int M){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a IssmDouble array");}
		void  SetValue(IssmDouble* pIssmDoublearray,int M,int N){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a IssmDouble array");}
		void  SetValue(int* intarray,int M){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a int array");}
		void  SetValue(int* pintarray,int M,int N){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a int array");}
		void  SetValue(Vector<IssmDouble>* vec){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a Vec");}
		void  SetValue(Matrix<IssmDouble>* mat){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a Mat");}
		void  SetValue(FILE* fid){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a FILE");}
		void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an array of matrices");}
		void  SetGradient(IssmDouble* poutput, int M, int N){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an IssmDouble");};
		/*}}}*/
};
#endif  /* _DOUBLEPARAM_H */
