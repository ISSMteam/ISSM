/*!\file KML_PolyStyle.cpp
 * \brief: implementation of the kml_polystyle object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Object.h"
#include "./KML_ColorStyle.h"
#include "./KML_PolyStyle.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_PolyStyle::KML_PolyStyle(){/*{{{*/

	fill      =true;
	outline   =true;

}
/*}}}*/
KML_PolyStyle::~KML_PolyStyle(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_PolyStyle::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_PolyStyle:\n");
	KML_ColorStyle::Echo();

	if(flag) _printf0_("          fill: " << fill << "\n");
	if(flag) _printf0_("       outline: " << outline << "\n");

	return;
}
/*}}}*/
void  KML_PolyStyle::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_PolyStyle::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_PolyStyle::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_PolyStyle:\n");
	KML_ColorStyle::DeepEcho(indent);

	if(flag) _printf0_(indent << "          fill: " << fill << "\n");
	if(flag) _printf0_(indent << "       outline: " << outline << "\n");

	return;
}
/*}}}*/
void  KML_PolyStyle::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<PolyStyle",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_ColorStyle::Write(filout,indent);

	fprintf(filout,"%s  <fill>%d</fill>\n",indent,fill);
	fprintf(filout,"%s  <outline>%d</outline>\n",indent,outline);

	fprintf(filout,"%s</PolyStyle>\n",indent);

	return;
}
/*}}}*/
void  KML_PolyStyle::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</PolyStyle",11)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_PolyStyle::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_PolyStyle::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<fill>"))
			KMLFileTokenParse(&fill      ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<outline>"))
			KMLFileTokenParse(&outline   ,
							  kstri,
							  fid);

		else if (!strncmp(kstri,"<",1))
			KML_ColorStyle::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
