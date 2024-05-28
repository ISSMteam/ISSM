/*! \file KML_LinearRing.h 
 *  \brief: header file for kml_linearring object
 */

#ifndef _KML_LINEARRING_H_
#define _KML_LINEARRING_H_

#define KML_LINEARRING_ALTMODE_LENGTH    18

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Geometry.h"
/*}}}*/

class KML_LinearRing: public KML_Geometry {

	public:

		bool     extrude;
		bool     tessellate;
		char     altmode[KML_LINEARRING_ALTMODE_LENGTH+1];
		int      ncoord;
		double  *coords;

		/*KML_LinearRing constructors, destructors {{{*/
		KML_LinearRing();
		~KML_LinearRing();
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
#endif  /* _KML_LINEARRING_H */
