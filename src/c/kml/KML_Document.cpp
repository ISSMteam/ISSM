/*!\file KML_Document.cpp
 * \brief: implementation of the kml_document object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Document.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Document::KML_Document(){/*{{{*/

	;

}
/*}}}*/
KML_Document::~KML_Document(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_Document::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Document:\n");
	KML_Container::Echo();

	return;
}
/*}}}*/
void  KML_Document::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Document::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Document::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Document:\n");
	KML_Container::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Document::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<Document",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Container::Write(filout,indent);

	fprintf(filout,"%s</Document>\n",indent);

	return;
}
/*}}}*/
void  KML_Document::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))) {
		if      (!strncmp(kstri,"</Document",10)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Document::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Document::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strncmp(kstri,"<",1))
			KML_Container::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for (ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
