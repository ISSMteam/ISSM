/*!\file KML_StyleSelector.cpp
 * \brief: implementation of the kml_styleselector abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Object.h"
#include "./KML_StyleSelector.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_StyleSelector::KML_StyleSelector(){/*{{{*/

	;

}
/*}}}*/
KML_StyleSelector::~KML_StyleSelector(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_StyleSelector::Echo(){/*{{{*/

	KML_Object::Echo();

	return;
}
/*}}}*/
void  KML_StyleSelector::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_StyleSelector::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_StyleSelector::DeepEcho(const char* indent){/*{{{*/

	KML_Object::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_StyleSelector::Write(FILE* filout,const char* indent){/*{{{*/

	KML_Object::Write(filout,indent);

	return;
}
/*}}}*/
void  KML_StyleSelector::Read(FILE* fid,char* kstr){/*{{{*/

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</StyleSelector",15))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_StyleSelector::Read -- Unexpected closing tag " << kstr << ".\n");}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_StyleSelector::Read -- Unexpected field \"" << kstr << "\".\n");}

	else if (!strncmp(kstr,"<",1))
		KML_Object::Read(fid,kstr);

	return;
}
/*}}}*/
