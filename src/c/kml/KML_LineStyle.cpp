/*!\file KML_LineStyle.cpp
 * \brief: implementation of the kml_linestyle object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_LineStyle.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_LineStyle::KML_LineStyle(){/*{{{*/

	width     =1.;

}
/*}}}*/
KML_LineStyle::~KML_LineStyle(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_LineStyle::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_LineStyle:\n");
	KML_ColorStyle::Echo();

	if(flag) _printf0_("         width: " << width << "\n");

	return;
}
/*}}}*/
void  KML_LineStyle::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_LineStyle::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_LineStyle::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_LineStyle:\n");
	KML_ColorStyle::DeepEcho(indent);

	if(flag) _printf0_(indent << "         width: " << width << "\n");

	return;
}
/*}}}*/
void  KML_LineStyle::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<LineStyle",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_ColorStyle::Write(filout,indent);

	fprintf(filout,"%s  <width>%g</width>\n",indent,width);

	fprintf(filout,"%s</LineStyle>\n",indent);

	return;
}
/*}}}*/
void  KML_LineStyle::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</LineStyle",11)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_LineStyle::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_LineStyle::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<width>"))
			KMLFileTokenParse(&width     ,
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
