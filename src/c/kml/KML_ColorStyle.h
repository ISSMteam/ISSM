/*! \file KML_ColorStyle.h 
 *  \brief: header file for kml_colorstyle abstract object
 */

#ifndef _KML_COLORSTYLE_H_
#define _KML_COLORSTYLE_H_

#define KML_COLORSTYLE_COLOR_LENGTH      8
#define KML_COLORSTYLE_COLORMODE_LENGTH  6

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_SubStyle.h"
/*}}}*/

class KML_ColorStyle: public KML_SubStyle {

	public:

		char  color[KML_COLORSTYLE_COLOR_LENGTH+1];
		char  colormode[KML_COLORSTYLE_COLORMODE_LENGTH+1];

		/*KML_ColorStyle constructors, destructors {{{*/
		KML_ColorStyle();
		~KML_ColorStyle();
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		void  Echo();
		void  DeepEcho();
		void  DeepEcho(const char* indent);
		void  Write(FILE* fid,const char* indent);
		void  Read(FILE* fid,char* kstr);
		int   Id(){_error_("Not implemented yet.");};
		void  Demarshall(char** pmarshalled_dataset){_error_("Not implemented yet.");};
		int   ObjectEnum(){_error_("Not implemented yet.");};
		Object* copy(){_error_("Not implemented yet.");};
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
		/*}}}*/

};
#endif  /* _KML_COLORSTYLE_H */
