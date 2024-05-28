/*! \file KML_Document.h 
 *  \brief: header file for kml_document object
 */

#ifndef _KML_DOCUMENT_H_
#define _KML_DOCUMENT_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Container.h"
class KML_Feature;
/*}}}*/

class KML_Document: public KML_Container {

	public:

		/*KML_Document constructors, destructors {{{*/
		KML_Document();
		~KML_Document();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
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
#endif  /* _KML_DOCUMENT_H */
