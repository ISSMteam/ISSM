/*! \file KML_Feature.h 
 *  \brief: header file for kml_feature abstract object
 */

#ifndef _KML_FEATURE_H_
#define _KML_FEATURE_H_

#define KML_FEATURE_NAME_LENGTH         80
#define KML_FEATURE_SNIPPET_LENGTH     160
#define KML_FEATURE_DESCRIPT_LENGTH   3200
#define KML_FEATURE_STYLEURL_LENGTH     80

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Object.h"
class KML_Style;
class DataSet;
/*}}}*/

class KML_Feature: public KML_Object {

	public:

		char  name[KML_FEATURE_NAME_LENGTH+1];
		bool  visibility;
		bool  open;
		char  snippet[KML_FEATURE_SNIPPET_LENGTH+1];
		char  descript[KML_FEATURE_DESCRIPT_LENGTH+1];
		char  styleurl[KML_FEATURE_STYLEURL_LENGTH+1];
		DataSet* style;

		/*KML_Feature constructors, destructors {{{*/
		KML_Feature();
		~KML_Feature();
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
#endif  /* _KML_FEATURE_H */
