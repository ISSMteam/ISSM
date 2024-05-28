#ifndef _INTARRAYINPUT_H_
#define _INTARRAYINPUT_H_

/*Headers:*/
#include "./Input.h"

class IntArrayInput: public Input{

	private:
		int         numberofelements_local;
		int*        N;
		int** values;

	public:
		/*IntArrayInput constructors, destructors: {{{*/
		IntArrayInput();
		IntArrayInput(int nbe_in);
		~IntArrayInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*IntArrayInput management:*/
		void SetInput(int row,int numinds,int* values_in);
		void GetArray(int row,int** pvalues,int* pN);
		void GetArrayPtr(int row,int** pvalues,int* pN);

};
#endif  /* _INTARRAYINPUT_H */
