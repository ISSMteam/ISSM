/*
 * GenericParam.h
 *
 *  Created on: Aug 29, 2012
 *      Author: utke
 */

#ifndef GENERICPARAM_H_
#define GENERICPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
        #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Param.h"
#include "../../shared/shared.h"

/**
 * here we have a class that holds an instance of P
 * but because it should live side by side with
 * the other instances derived from Param it - unfortunately -
 * inherits all the accessors that are useless in this context
 */
template <class P> class GenericParam: public Param{

        private:
                P myP;
                int myEnumVal;

        public:
                /*GenericParam constructors, destructors: {{{*/
                GenericParam(int enumVal) : myEnumVal(enumVal){};
                GenericParam(P Pin, int enumVal) : myP(Pin),myEnumVal(enumVal){};
                ~GenericParam(){};
                /*}}}*/
                /*Object virtual functions definitions:{{{ */
                // unfortunately,  having to implement such a printer method implies
                // that any structured P must provide the friend << operator
                Param* copy() { return new GenericParam<P>(*this); };
                void  DeepEcho() {
                  _printf_("GenericParam:\n");
                  _printf_("   enum:  " << myEnumVal << " (" << EnumToStringx(myEnumVal) << ")\n");
                  _printf_("   value: " << myP << "\n");;
                }
                void  Echo() {DeepEcho();};
                int   Id(){ return -1; };

                // the "copy"  has to implement the base class abstract function
                // but I would prefer to drop this not to hide a "new" in here because
                // it does not clarify  ownership of the newed up instance...
                // use the default copy constructor instead
                void Marshall(MarshallHandle* marshallhandle){
						 if(this->InstanceEnum()!=FemModelCommEnum){
							 _printf_("   WARNING: parameter "<<EnumToStringx(this->myEnumVal)<<" is a GenericParam and cannot be marshalled\n");
						 }
						 /*Nothing for now*/
					 }
                int   ObjectEnum() {return GenericParamEnum;};

                /*}}}*/
                /*Param virtual function definitions: {{{*/
                P& GetParameterValue() { return myP;}
                const P& GetParameterValue()const { return myP;};
                int   InstanceEnum(){return myEnumVal;}

                // none of these apply ...
                void  GetParameterValue(bool* pbool){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a bool");}
                void  GetParameterValue(int* pinteger){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return an integer");}
                void  GetParameterValue(int** pintarray,int* pM){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return an array of integers");}
                void  GetParameterValue(int** pintarray,int* pM,int* pN){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return an array of integers");}
                void  GetParameterValue(IssmDouble* pIssmDouble){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble");}
                void  GetParameterValue(IssmDouble* pdouble,IssmDouble time){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble for a given time");}
                void  GetParameterValue(IssmDouble* pdouble,IssmDouble time, int timestepping, IssmDouble dt){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble for a given time");}
                void  GetParameterValue(char** pstring){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a string");}
                void  GetParameterValue(char*** pstringarray,int* pM){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a string array");}
                void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble array");}
                void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, int* pN){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble array");}
                void  GetParameterValue(IssmDouble** pIssmDoublearray,int* pM, const char* data){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a IssmDouble array");}
                void  GetParameterValue(IssmDouble*** parray, int* pM,int** pmdims, int** pndims){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a matrix array");}
                void  GetParameterValue(Vector<IssmDouble>** pvec){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a Vec");}
                void  GetParameterValue(Matrix<IssmDouble>** pmat){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a Mat");}
                void  GetParameterValue(FILE** pfid){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a FILE");}
                void  GetParameterValue(DataSet** pdataset){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot return a DataSet");}

                void  SetEnum(int enum_in){this->myEnumVal = enum_in;};
                void  SetValue(bool boolean){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a bool");}
                void  SetValue(int integer){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an integer");}
                void  SetValue(int* intarray,int M){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an int array");}
                void  SetValue(int* intarray,int M,int N){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an int array");}
                void  SetValue(IssmDouble scalar){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an IssmDouble");}
                void  SetValue(char* string){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a string");}
                void  SetValue(char** stringarray,int M){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a string array");}
                void  SetValue(IssmDouble* IssmDoublearray){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a IssmDouble array");}
                void  SetValue(IssmDouble* IssmDoublearray,int M){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a IssmDouble array");}
                void  SetValue(IssmDouble* pIssmDoublearray,int M,int N){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a IssmDouble array");}
                void  SetValue(Vector<IssmDouble>* vec){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a Vec");}
                void  SetValue(Matrix<IssmDouble>* mat){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a Mat");}
                void  SetValue(FILE* fid){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold a FILE");}
                void  SetValue(IssmDouble** array, int M, int* mdim_array, int* ndim_array){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an array of matrices");}
                void  SetGradient(IssmDouble* poutput, int M, int N){_error_("Param "<< EnumToStringx(myEnumVal) << " cannot hold an IssmDouble");};

                /*}}}*/
};

#endif /* GENERICPARAM_H_ */
