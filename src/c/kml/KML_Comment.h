/*! \file KML_Comment.h 
 *  \brief: header file for kml_comment object
 */

#ifndef _KML_COMMENT_H_
#define _KML_COMMENT_H_

/*Headers:{{{*/
#include "../shared/shared.h"
#include "../datastructures/datastructures.h"
class DataSet;
/*}}}*/

class KML_Comment: public Object {

	public:

		char* name;
		char* value;

		/*KML_Comment constructors, destructors {{{*/
		KML_Comment();
		~KML_Comment();
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		virtual void  Echo();
		virtual void  DeepEcho();
		virtual void  DeepEcho(const char* indent);
		int   Id(){_error_("Not implemented yet.");};
		int   ObjectEnum(){_error_("Not implemented yet.");};
		Object* copy(){_error_("Not implemented yet.");};
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
		/*}}}*/

		/*virtual functions: */
		void  Write(FILE* fid,const char* indent);
		void  Read(FILE* fid,char* kstr);
		void  Alloc(const char* valuei);
		void  Add(DataSet* commnt);
		void  Get(char** pvalueo);

};
#endif  /* _KML_COMMENT_H */
