#ifndef _DOUBLEINPUT2_H_
#define _DOUBLEINPUT2_H_

/*Headers:*/
#include "./Input.h"

class DoubleInput: public Input{

	private:
		int   size;
		IssmDouble*  values;

	public:
		/*DoubleInput constructors, destructors: {{{*/
		DoubleInput();
		DoubleInput(int size_in);
		~DoubleInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*DoubleInput management: {{{*/
		void GetInput(IssmDouble* pvalue,int index);
		void SetInput(int index,IssmDouble value);
		/*}}}*/

};
#endif  /* _BOOLINPUT_H */
