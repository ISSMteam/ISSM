/*! \file KML_Icon.h 
 *  \brief: header file for kml_icon object
 */

#ifndef _KML_ICON_H_
#define _KML_ICON_H_

#define KML_ICON_HREF_LENGTH      800
#define KML_ICON_REFMODE_LENGTH    10
#define KML_ICON_VREFMODE_LENGTH    9
#define KML_ICON_VFORMAT_LENGTH   800
#define KML_ICON_HQUERY_LENGTH    800

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
/*}}}*/

class KML_Icon: public KML_Object {

	public:

		char  href[KML_ICON_HREF_LENGTH+1];
		char  refmode[KML_ICON_REFMODE_LENGTH+1];
		float refint;
		char  vrefmode[KML_ICON_VREFMODE_LENGTH+1];
		float vreftime;
		float vboundsc;
		char  vformat[KML_ICON_VFORMAT_LENGTH+1];
		char  hquery[KML_ICON_HQUERY_LENGTH+1];

		/*KML_Icon constructors, destructors {{{*/
		KML_Icon();
		~KML_Icon();
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
#endif  /* _KML_ICON_H */
