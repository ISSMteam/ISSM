/*! \file KML_StyleSelector.h 
 *  \brief: header file for kml_styleselector abstract object
 */

#ifndef _KML_STYLESELECTOR_H_
#define _KML_STYLESELECTOR_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
/*}}}*/

class KML_StyleSelector: public KML_Object {

	public:

		/*KML_StyleSelector constructors, destructors {{{*/
		KML_StyleSelector();
		~KML_StyleSelector();
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
#endif  /* _KML_STYLESELECTOR_H */
