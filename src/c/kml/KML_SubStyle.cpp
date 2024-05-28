/*!\file KML_SubStyle.cpp
 * \brief: implementation of the kml_substyle abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_SubStyle.h"
#include "./KML_Object.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_SubStyle::KML_SubStyle(){/*{{{*/

	;

}
/*}}}*/
KML_SubStyle::~KML_SubStyle(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_SubStyle::Echo(){/*{{{*/

	KML_Object::Echo();

	return;
}
/*}}}*/
void  KML_SubStyle::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_SubStyle::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_SubStyle::DeepEcho(const char* indent){/*{{{*/

	KML_Object::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_SubStyle::Write(FILE* filout,const char* indent){/*{{{*/

	KML_Object::Write(filout,indent);

	return;
}
/*}}}*/
void  KML_SubStyle::Read(FILE* fid,char* kstr){/*{{{*/

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</SubStyle",10))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_SubStyle::Read -- Unexpected closing tag " << kstr << ".\n");}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_SubStyle::Read -- Unexpected field \"" << kstr << "\".\n");}

	else if (!strncmp(kstr,"<",1))
		KML_Object::Read(fid,kstr);

	return;
}
/*}}}*/
