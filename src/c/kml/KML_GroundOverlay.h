/*! \file KML_GroundOverlay.h 
 *  \brief: header file for kml_groundoverlay object
 */

#ifndef _KML_GROUNDOVERLAY_H_
#define _KML_GROUNDOVERLAY_H_

#define KML_GROUNDOVERLAY_ALTMODE_LENGTH    18

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Overlay.h"
class KML_LatLonBox;
/*}}}*/

class KML_GroundOverlay: public KML_Overlay {

	public:

		double altitude;
		char  altmode[KML_GROUNDOVERLAY_ALTMODE_LENGTH+1];
		KML_LatLonBox* llbox;

		/*KML_GroundOverlay constructors, destructors {{{*/
		KML_GroundOverlay();
		~KML_GroundOverlay();
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
#endif  /* _KML_GROUNDOVERLAY_H */
