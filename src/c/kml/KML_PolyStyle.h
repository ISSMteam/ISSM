/*! \file KML_PolyStyle.h 
 *  \brief: header file for kml_polystyle object
 */

#ifndef _KML_POLYSTYLE_H_
#define _KML_POLYSTYLE_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_ColorStyle.h"
/*}}}*/

class KML_PolyStyle: public KML_ColorStyle {

	public:

		int   fill;
		int   outline;

		/*KML_PolyStyle constructors, destructors {{{*/
		KML_PolyStyle();
		~KML_PolyStyle();
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
#endif  /* _KML_POLYSTYLE_H */
