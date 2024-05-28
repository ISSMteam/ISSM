/*! \file KML_Point.h 
 *  \brief: header file for kml_point object
 */

#ifndef _KML_POINT_H_
#define _KML_POINT_H_

#define KML_POINT_ALTMODE_LENGTH    18

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Geometry.h"
/*}}}*/

class KML_Point: public KML_Geometry {

	public:

		bool  extrude;
		char  altmode[KML_POINT_ALTMODE_LENGTH+1];
		double coords[3];

		/*KML_Point constructors, destructors {{{*/
		KML_Point();
		~KML_Point();
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		void  Echo();
		void  DeepEcho();
		void  DeepEcho(const char* indent);
		void  Write(FILE* fid,const char* indent);
		void  Read(FILE* fid,char* kstr);
		void  WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp);
		int   Id(){_error_("Not implemented yet.");};
		int   ObjectEnum(){_error_("Not implemented yet.");};
		Object* copy(){_error_("Not implemented yet.");};
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
		/*}}}*/

};
#endif  /* _KML_POINT_H */
