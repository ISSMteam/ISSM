/*! \file KML_LatLonBox.h 
 *  \brief: header file for kml_latlonbox object
 */

#ifndef _KML_LATLONBOX_H_
#define _KML_LATLONBOX_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
/*}}}*/

class KML_LatLonBox: public KML_Object {

	public:

		double north;
		double south;
		double east;
		double west;
		double rotation;

		/*KML_LatLonBox constructors, destructors {{{*/
		KML_LatLonBox();
		~KML_LatLonBox();
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
#endif  /* _KML_LATLONBOX_H */
