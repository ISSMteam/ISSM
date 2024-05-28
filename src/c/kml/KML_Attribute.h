/*! \file KML_Attribute.h 
 *  \brief: header file for kml_attribute object
 */

#ifndef _KML_ATTRIBUTE_H_
#define _KML_ATTRIBUTE_H_

/*Headers:{{{*/
#include "../shared/shared.h"
#include "../datastructures/datastructures.h"
/*}}}*/

class KML_Attribute: public Object {

	public:

		char* name;
		char* value;

		/*KML_Attribute constructors, destructors {{{*/
		KML_Attribute();
		~KML_Attribute();
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		virtual void  Echo();
		virtual void  DeepEcho();
		virtual void  DeepEcho(const char* indent);
		int   Id(){_error_("Not implemented yet.");};
		int   ObjectEnum(){_error_("Not implemented yet.");};
		Object* copy(){_error_("Not implemented yet.");};
		void Marshall(MarshallHandle* marshallhandle);
		/*}}}*/

		/*virtual functions: */
		void  Write(FILE* fid,const char* indent);
		void  Read(FILE* fid,char* kstr);
		void  Alloc(const char* namei,const char* valuei);
		void  Add(DataSet* attrib);
		void  Get(char** pvalueo,char* deflt);

};
#endif  /* _KML_ATTRIBUTE_H */
