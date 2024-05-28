/*! \file KML_Polygon.h 
 *  \brief: header file for kml_polygon object
 */

#ifndef _KML_POLYGON_H_
#define _KML_POLYGON_H_

#define KML_POLYGON_ALTMODE_LENGTH    18

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Geometry.h"
class KML_LinearRing;
class DataSet;
/*}}}*/

class KML_Polygon: public KML_Geometry {

	public:

		bool  extrude;
		bool  tessellate;
		char  altmode[KML_POLYGON_ALTMODE_LENGTH+1];
		DataSet* outer;
		DataSet* inner;

		/*KML_Polygon constructors, destructors {{{*/
		KML_Polygon();
		~KML_Polygon();
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
#endif  /* _KML_POLYGON_H */
