/*! \file KML_Object.h 
 *  \brief: header file for kml_object abstract object
 */

#ifndef _KML_OBJECT_H_
#define _KML_OBJECT_H_

/*Headers:{{{*/
#include "../shared/shared.h"
#include "../datastructures/datastructures.h"
/*}}}*/

class KML_Object: public Object {

	public:

		DataSet* attrib;
		DataSet* commnt;
		DataSet* kmlobj;

		/*KML_Object constructors, destructors {{{*/
		KML_Object();
		~KML_Object();
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
		virtual void  Write(FILE* fid,const char* indent)=0;
		virtual void  Read(FILE* fid,char* kstr)=0;
		virtual void  WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp);
		virtual void  AddAttrib(const char* name,const char* value);
		virtual void  WriteAttrib(FILE* fid,const char* indent);
		virtual void  AddCommnt(int ncom,char** pcom);
		virtual void  AddCommnt(char* value);
		virtual void  WriteCommnt(FILE* fid,const char* indent);

};
#endif  /* _KML_OBJECT_H */
