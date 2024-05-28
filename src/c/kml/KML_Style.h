/*! \file KML_Style.h 
 *  \brief: header file for kml_style object
 */

#ifndef _KML_STYLE_H_
#define _KML_STYLE_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_StyleSelector.h"
class KML_LineStyle;
class KML_PolyStyle;
/*}}}*/

class KML_Style: public KML_StyleSelector {

	public:

		void* icon;
		void* label;
		KML_LineStyle* line;
		KML_PolyStyle* poly;
		void* balloon;
		void* list;

		/*KML_Style constructors, destructors {{{*/
		KML_Style();
		~KML_Style();
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
#endif  /* _KML_STYLE_H */
