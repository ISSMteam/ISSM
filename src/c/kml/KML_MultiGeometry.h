/*! \file KML_MultiGeometry.h 
 *  \brief: header file for kml_multigeometry object
 */

#ifndef _KML_MULTIGEOMETRY_H_
#define _KML_MULTIGEOMETRY_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Geometry.h"
class KML_Geometry;
class DataSet;
/*}}}*/

class KML_MultiGeometry: public KML_Geometry {

	public:

		DataSet* geometry;

		/*KML_MultiGeometry constructors, destructors {{{*/
		KML_MultiGeometry();
		~KML_MultiGeometry();
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
#endif  /* _KML_MULTIGEOMETRY_H */
