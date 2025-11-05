/*! \file VectorParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _VECTORPARAM_H_
#define _VECTORPARAM_H_

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

class VectorParam: public Param{

	private: 
		/*just hold 3 values for 3 vertices: */
		int enum_type;
		Vector<IssmDouble>* value;

	public:
		/*VectorParam constructors, destructors: {{{*/
		VectorParam();
		VectorParam(int enum_type,Vector<IssmDouble>* value);
		~VectorParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!"); };
		int   ObjectEnum();
		/*}}}*/
		/*Param vritual function definitions: {{{*/
		void  GetParameterValue(bool* pbool){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a bool");}
		void  GetParameterValue(int* pinteger){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an integer");}
		void  GetParameterValue(int** pintarray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an array of integers");}
		void  GetParameterValue(int** pintarray,int* pM,int* pN){_error_("Param "<< EnumToStringx(enum_type) << " cannot return an array of integers");}
		void  GetParameterValue(IssmDouble* pIssmDouble){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble");}
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble for a given time");}
		void  GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble for a given time");}
		void  GetParameterValue(char** pstring){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string");}
		void  GetParameterValue(char*** pstringarray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a string array");}
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble array");}
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, int* pN){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble array");}
		void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a IssmDouble array");}
		void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a matrix array");}
		void  GetParameterValue(Matrix<IssmDouble>** pmat){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a Mat");}
		void  GetParameterValue(Vector<IssmDouble>** poutput);
		void  GetParameterValue(FILE** pfid){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot return a FILE");}
		void  GetParameterValue(DataSet** pdataset){_error_("Param "<< EnumToStringx(enum_type) << " cannot return a DataSet");}
		int   InstanceEnum(){return enum_type;}

		void  SetEnum(int enum_in){this->enum_type = enum_in;};
		void  SetValue(bool boolean){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a boolean");}
		void  SetValue(int integer){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold an integer");}
		void  SetValue(IssmDouble scalar){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a scalar");}
		void  SetValue(char* string){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a string");}
		void  SetValue(char** stringarray,int M){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a string array");}
		void  SetValue(IssmDouble* IssmDoublearray){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold a IssmDouble array");}
		void  SetValue(IssmDouble* IssmDoublearray,int M){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a IssmDouble array");}
		void  SetValue(IssmDouble* pIssmDoublearray,int M,int N){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a IssmDouble array");}
		void  SetValue(int* intarray,int M){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a int array");}
		void  SetValue(int* pintarray,int M,int N){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a int array");}
		void  SetValue(Vector<IssmDouble>* vec);
		void  SetValue(Matrix<IssmDouble>* mat){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a Mat");}
		void  SetValue(FILE* fid){_error_("Vector of enum " << enum_type << " (" << EnumToStringx(enum_type) << ") cannot hold a FILE");}
		void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an array of matrices");}
		void  SetGradient(IssmDouble* poutput, int M, int N){_error_("Param "<< EnumToStringx(enum_type) << " cannot hold an IssmDouble");};
		/*}}}*/
};
#endif  /* _VECTORPARAM_H */
