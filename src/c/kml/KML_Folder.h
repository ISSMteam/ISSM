/*! \file KML_Folder.h 
 *  \brief: header file for kml_folder object
 */

#ifndef _KML_FOLDER_H_
#define _KML_FOLDER_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Container.h"
class KML_Feature;
/*}}}*/

class KML_Folder: public KML_Container {

	public:

		/*KML_Folder constructors, destructors {{{*/
		KML_Folder();
		~KML_Folder();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
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
#endif  /* _KML_FOLDER_H */
