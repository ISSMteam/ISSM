#ifndef _INTINPUT2_H_
#define _INTINPUT2_H_

/*Headers:*/
#include "./Input.h"

class IntInput: public Input{

	private:
		int   size;
		int*  values;

	public:
		/*IntInput constructors, destructors: {{{*/
		IntInput();
		IntInput(int size_in);
		~IntInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*IntInput management: {{{*/
		void GetInput(int* pvalue,int index);
		void SetInput(int index,int value);
		/*}}}*/

};
#endif  /* _BOOLINPUT_H */
