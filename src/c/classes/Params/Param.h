/*!\file:  Param.h
 * \brief abstract class for Param object
 */ 

#ifndef _PARAM_H_
#define _PARAM_H_

/*Headers:*/
/*{{{*/

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../datastructures/datastructures.h"
#include "../Node.h"
/*}}}*/

class Param: public Object{

	protected:
		int enum_type;

	public: 

		/*Non virtual functions*/
		int   Id(){return -1;}
		int   InstanceEnum(){return enum_type;}
		void  SetEnum(int enum_in){this->enum_type = enum_in;}

		/*Virtual functions:*/
		virtual        ~Param(){};
		virtual void   DeepEcho()=0;
		virtual Param* copy()=0;
		virtual void   Echo()=0;
		virtual void  GetParameterValue(bool* pbool){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a bool");}
		virtual void  GetParameterValue(int* pinteger){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an integer");}
		virtual void  GetParameterValue(int** pintarray,int* pM){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an array of integers");}
		virtual void  GetParameterValue(int** pintarray,int* pM,int* pN){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an array of integers");}
		virtual void  GetParameterValue(IssmDouble* pIssmDouble){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble");}
		virtual void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble for a given time");}
		virtual void  GetParameterValue(IssmDouble* pdouble,IssmDouble time,int timestepping,IssmDouble dt){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble for a given time");}
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, int column, IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, IssmDouble time, int timestepping, IssmDouble dt){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, int column, IssmDouble time, int timestepping, IssmDouble dt){_error_("not implemented yet");};
		virtual void  GetParameterValue(char** pstring){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a string");}
		virtual void  GetParameterValue(char*** pstringarray,int* pM){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a string array");}
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble array");}
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble starttime,IssmDouble endtime){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble array");}
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return an IssmDouble array");}
		virtual void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a matrix array");}
		virtual void  GetParameterValue(Vector<IssmDouble>** pvec){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a Vec");}
		virtual void  GetParameterValue(Matrix<IssmDouble>** pmat){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a Mat");}
		virtual void  GetParameterValue(FILE** pfid){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a FILE");}
		virtual void  GetParameterValue(DataSet** pdataset){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot return a DataSet");}
		virtual void  Marshall(MarshallHandle* marshallhandle)=0;
		virtual int   ObjectEnum()=0;
		virtual void  SetValue(bool boolean){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a bool");}
		virtual void  SetValue(int integer){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an int");}
		virtual void  SetValue(IssmDouble scalar){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an IssmDouble");}
		virtual void  SetValue(char* string){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a string");}
		virtual void  SetValue(char** stringarray,int M){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a string array");}
		virtual void  SetValue(DataSet* dataset){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a DataSet");}
		virtual void  SetValue(IssmDouble* IssmDoublearray){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an IssmDouble array");}
		virtual void  SetValue(IssmDouble* IssmDoublearray,int M){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an IssmDouble array");}
		virtual void  SetValue(IssmDouble* pIssmDoublearray,int M,int N){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an IssmDouble array");}
		virtual void  SetValue(int* intarray,int M){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an int array");}
		virtual void  SetValue(int* pintarray,int M,int N){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an int array");}
		virtual void  SetValue(Vector<IssmDouble>* vec){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a Vec");}
		virtual void  SetValue(Matrix<IssmDouble>* mat){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a Mat");}
		virtual void  SetValue(FILE* fid){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold a FILE");}
		virtual void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array){_error_("Param "<< EnumToStringx(this->enum_type) << " cannot hold an array of matrices");}
		virtual void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int N,const char* data,int offset){_error_("not implemented yet");};
};
#endif
