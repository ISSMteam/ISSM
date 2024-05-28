/*! \file KML_File.h 
 *  \brief: header file for kml_file object
 */

#ifndef _KML_FILE_H_
#define _KML_FILE_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"

#include "./KML_Feature.h"
class DataSet;
/*}}}*/

class KML_File: public KML_Object {

	public:

		/*KML_File constructors, destructors {{{*/
		KML_File();
		~KML_File();
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
		void Marshall(MarshallHandle* marshallhandle);
		/*}}}*/

};
#endif  /* _KML_FILE_H */
