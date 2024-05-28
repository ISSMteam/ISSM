/*! \file KML_SubStyle.h 
 *  \brief: header file for kml_substyle abstract object
 */

#ifndef _KML_SUBSTYLE_H_
#define _KML_SUBSTYLE_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
/*}}}*/

class KML_SubStyle: public KML_Object {

	public:

		/*KML_SubStyle constructors, destructors {{{*/
		KML_SubStyle();
		~KML_SubStyle();
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		void  Echo();
		void  DeepEcho();
		void  DeepEcho(const char* indent);
		void  Write(FILE* fid,const char* indent);
		void  Read(FILE* fid,char* kstr);
		int   Id(){_error_("Not implemented yet.");};
		int   ObjectEnum(){_error_("Not implemented yet.");};
		Object* copy(){_error_("Not implemented yet.");};
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
		/*}}}*/

};
#endif  /* _KML_SUBSTYLE_H */
