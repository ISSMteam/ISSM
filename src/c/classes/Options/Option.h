/*! \file Option.h 
 *  \brief: header file for option abstract object
 */

#ifndef _OPTIONOBJECT_H_
#define _OPTIONOBJECT_H_

/*Headers:{{{*/
#include "../../shared/shared.h"
#include "../../datastructures/datastructures.h"
/*}}}*/

class Option: public Object {

	public:

		/*Option constructors, destructors*/
		Option(){};
		~Option(){};

		/*Object virtual functions definitions*/
		Object       *copy(){_error_("Not implemented yet"); };
		virtual void  DeepEcho()= 0;
		virtual void  DeepEcho(char  *indent)=0;
		virtual void  Echo()= 0;
		int           Id(){_error_("Not implemented yet"); };
		void          Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!"); };
		int           ObjectEnum(){return OptionEnum;};

		/*virtual functions: */
		virtual char* Name()=0;

};
#endif  /* _OPTIONOBJECT_H */
