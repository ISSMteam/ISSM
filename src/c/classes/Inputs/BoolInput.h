#ifndef _BOOLINPUT2_H_
#define _BOOLINPUT2_H_

/*Headers:*/
#include "./Input.h"

class BoolInput: public Input{

	private:
		int   size;
		bool* values;

	public:
		/*BoolInput constructors, destructors: {{{*/
		BoolInput();
		BoolInput(int size_in);
		~BoolInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*BoolInput management: {{{*/
		void GetInput(bool* pvalue,int index);
		void SetInput(int index,bool value);
		/*}}}*/
		/*numerics: {{{*/
		/*}}}*/

};
#endif  /* _BOOLINPUT_H */
