/*! \file KML_Container.h 
 *  \brief: header file for kml_container abstract object
 */

#ifndef _KML_CONTAINER_H_
#define _KML_CONTAINER_H_

/*Headers:*/
/*{{{*/
#include "../shared/shared.h"
#include "./KML_Feature.h"
class DataSet;
/*}}}*/

class KML_Container: public KML_Feature {

	public:

		DataSet* feature;

		/*KML_Container constructors, destructors {{{*/
		KML_Container();
		~KML_Container();
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
#endif  /* _KML_CONTAINER_H */
