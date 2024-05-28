/*!\file KML_Geometry.cpp
 * \brief: implementation of the kml_geometry abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Geometry.h"
#include "./KML_Object.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Geometry::KML_Geometry(){/*{{{*/

	;

}
/*}}}*/
KML_Geometry::~KML_Geometry(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_Geometry::Echo(){/*{{{*/

	this->KML_Object::Echo();

	return;
}
/*}}}*/
void  KML_Geometry::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Geometry::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Geometry::DeepEcho(const char* indent){/*{{{*/

	this->KML_Object::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Geometry::Write(FILE* filout,const char* indent){/*{{{*/

	KML_Object::Write(filout,indent);

	return;
}
/*}}}*/
void  KML_Geometry::Read(FILE* fid,char* kstr){/*{{{*/

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</Geometry",10))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_Geometry::Read -- Unexpected closing tag " << kstr << ".\n");}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_Geometry::Read -- Unexpected field \"" << kstr << "\".\n");}

	else if (!strncmp(kstr,"<",1))
		KML_Object::Read(fid,kstr);

	return;
}
/*}}}*/
