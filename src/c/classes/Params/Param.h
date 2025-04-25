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

	public: 
		virtual        ~Param(){};

		/*Virtual functions:*/
		virtual void  DeepEcho()=0;
		virtual Param* copy()=0;
		virtual void  Echo()=0;
		virtual void  GetParameterValue(bool* pbool)=0;
		virtual void  GetParameterValue(int* pinteger)=0;
		virtual void  GetParameterValue(int** pintarray,int* pM)=0;
		virtual void  GetParameterValue(int** pintarray,int* pM,int* pN)=0;
		virtual void  GetParameterValue(IssmDouble* pIssmDouble)=0;
		virtual void  GetParameterValue(IssmDouble* pdouble,IssmDouble time)=0;
		virtual void  GetParameterValue(IssmDouble* pdouble,IssmDouble time,int timestepping,IssmDouble dt)=0;
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, int column, IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, IssmDouble time, int timestepping, IssmDouble dt){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble* pdouble,int row, int column, IssmDouble time, int timestepping, IssmDouble dt){_error_("not implemented yet");};
		virtual void  GetParameterValue(char** pstring)=0;
		virtual void  GetParameterValue(char*** pstringarray,int* pM)=0;
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM)=0;
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble time){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble starttime,IssmDouble endtime){_error_("not implemented yet");};
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM,int* pN)=0;
		virtual void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data)=0;
		virtual void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims)=0;
		virtual void  GetParameterValue(Vector<IssmDouble>** pvec)=0;
		virtual void  GetParameterValue(Matrix<IssmDouble>** pmat)=0;
		virtual void  GetParameterValue(FILE** pfid)=0;
		virtual void  GetParameterValue(DataSet** pdataset)=0;
		virtual int   InstanceEnum()=0;
		virtual void  Marshall(MarshallHandle* marshallhandle)=0;
		virtual int   ObjectEnum()=0;

		virtual void  SetEnum(int enum_in)=0;
		virtual void  SetValue(bool boolean)=0;
		virtual void  SetValue(int integer)=0;
		virtual void  SetValue(IssmDouble scalar)=0;
		virtual void  SetValue(char* string)=0;
		virtual void  SetValue(char** stringarray,int M)=0;
		virtual void  SetValue(DataSet* dataset){_error_("not implemented yet");};
		virtual void  SetValue(IssmDouble* IssmDoublearray,int M)=0;
		virtual void  SetValue(IssmDouble* IssmDoublearray)=0;
		virtual void  SetValue(IssmDouble* pIssmDoublearray,int M,int N)=0;
		virtual void  SetValue(int* intarray,int M)=0;
		virtual void  SetValue(int* pintarray,int M,int N)=0;
		virtual void  SetValue(Vector<IssmDouble>* vec)=0;
		virtual void  SetValue(Matrix<IssmDouble>* mat)=0;
		virtual void  SetValue(FILE* fid)=0;
		virtual void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array)=0;
		virtual void  SetGradient(IssmDouble* poutput, int M, int N)=0;
		virtual void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_index,int N,const char* data,int offset){_error_("not implemented yet");};
};
#endif
