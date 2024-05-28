/*! \file KML_Unknown.h 
 *  \brief: header file for kml_unknown object
 */

#ifndef _KML_UNKNOWN_H_
#define _KML_UNKNOWN_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
/*}}}*/

class KML_Unknown: public KML_Object {

	public:

		char* name;
		char* value;

		/*KML_Unknown constructors, destructors {{{*/
		KML_Unknown();
		~KML_Unknown();
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
#endif  /* _KML_UNKNOWN_H */
