/*! \file KML_LineStyle.h 
 *  \brief: header file for kml_linestyle object
 */

#ifndef _KML_LINESTYLE_H_
#define _KML_LINESTYLE_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_ColorStyle.h"
/*}}}*/

class KML_LineStyle: public KML_ColorStyle {

	public:

		float width;

		/*KML_LineStyle constructors, destructors {{{*/
		KML_LineStyle();
		~KML_LineStyle();
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
#endif  /* _KML_LINESTYLE_H */
