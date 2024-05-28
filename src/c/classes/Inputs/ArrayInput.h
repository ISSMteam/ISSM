#ifndef _ARRAYINPUT2_H_
#define _ARRAYINPUT2_H_

/*Headers:*/
#include "./Input.h"

class ArrayInput: public Input{

	private:
		int         numberofelements_local;
		int*        N;
		IssmDouble** values;

	public:
		/*ArrayInput constructors, destructors: {{{*/
		ArrayInput();
		ArrayInput(int nbe_in);
		~ArrayInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*ArrayInput management:*/
		void SetInput(int row,int numinds,IssmDouble* values_in);
		void GetArray(int row,IssmDouble** pvalues,int* pN);
		void GetArrayPtr(int row,IssmDouble** pvalues,int* pN);

};
#endif  /* _ARRAYINPUT_H */
