#ifndef _INTARRAYINPUT_H_
#define _INTARRAYINPUT_H_

/*Headers:*/
#include <cstdint>
#include "./Input.h"

template<typename T>
class TIntArrayInput: public Input{

	private:
		int         numberofelements_local;
		int*        N;
		T**         values;

	public:
		/*TIntArrayInput constructors, destructors:*/
		TIntArrayInput();
		TIntArrayInput(int nbe_in);
		~TIntArrayInput();

		/*Object virtual functions definitions:*/
		Input *copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum();

		/*TIntArrayInput management:*/
		void SetInput(int row,int numinds,T* values_in);
		void GetArray(int row,T** pvalues,int* pN);
		void GetArrayPtr(int row,T** pvalues,int* pN);

};

/*ObjectEnum is specialized for each supported type (see IntArrayInput.cpp)*/
template<> int TIntArrayInput<int>::ObjectEnum();
template<> int TIntArrayInput<int16_t>::ObjectEnum();

/*Default (historical) 32-bit version*/
typedef TIntArrayInput<int> IntArrayInput;

/*Compact 16-bit version (e.g. for quantized indices)*/
typedef TIntArrayInput<int16_t> Int16ArrayInput;

#endif  /* _INTARRAYINPUT_H */
