/*! \file KML_Overlay.h 
 *  \brief: header file for kml_overlay abstract object
 */

#ifndef _KML_OVERLAY_H_
#define _KML_OVERLAY_H_

#define KML_OVERLAY_COLOR_LENGTH  8

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Feature.h"
class KML_Icon;
/*}}}*/

class KML_Overlay: public KML_Feature {

	public:

		char  color[KML_OVERLAY_COLOR_LENGTH+1];
		int   draword;
		KML_Icon* icon;

		/*KML_Overlay constructors, destructors {{{*/
		KML_Overlay();
		~KML_Overlay();
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
#endif  /* _KML_OVERLAY_H */
