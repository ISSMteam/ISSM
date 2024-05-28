/*!\file KML_LatLonBox.cpp
 * \brief: implementation of the kml_feature abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_LatLonBox.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_LatLonBox::KML_LatLonBox(){/*{{{*/

	north     = 0.;
	south     = 0.;
	east      = 0.;
	west      = 0.;
	rotation  = 0.;

}
/*}}}*/
KML_LatLonBox::~KML_LatLonBox(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_LatLonBox::Echo(){/*{{{*/

	_printf_("KML_LatLonBox:\n");
	KML_Object::Echo();

	_printf_("         north: " << north << "\n");
	_printf_("         south: " << south << "\n");
	_printf_("          east: " << east << "\n");
	_printf_("          west: " << west << "\n");
	_printf_("      rotation: " << rotation << "\n");
}
/*}}}*/
void  KML_LatLonBox::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_LatLonBox::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_LatLonBox::DeepEcho(const char* indent){/*{{{*/

	_printf_(indent << "KML_LatLonBox:\n");
	KML_Object::DeepEcho(indent);

	_printf_("         north: " << north << "\n");
	_printf_("         south: " << south << "\n");
	_printf_("          east: " << east << "\n");
	_printf_("          west: " << west << "\n");
	_printf_("      rotation: " << rotation << "\n");
}
/*}}}*/
void  KML_LatLonBox::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<LatLonBox",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Object::Write(filout,indent);

	fprintf(filout,"%s  <north>%0.16g</north>\n",indent,north);
	fprintf(filout,"%s  <south>%0.16g</south>\n",indent,south);
	fprintf(filout,"%s  <east>%0.16g</east>\n",indent,east);
	fprintf(filout,"%s  <west>%0.16g</west>\n",indent,west);
	fprintf(filout,"%s  <rotation>%0.16g</rotation>\n",indent,rotation);

	fprintf(filout,"%s</LatLonBox>\n",indent);

	return;
}
/*}}}*/
void  KML_LatLonBox::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</LatLonBox",11)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_LatLonBox::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_LatLonBox::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<north>"))
			KMLFileTokenParse(&north     ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<south>"))
			KMLFileTokenParse(&south     ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<east>"))
			KMLFileTokenParse(&east      ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<west>"))
			KMLFileTokenParse(&west      ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<rotation>"))
			KMLFileTokenParse(&rotation  ,
							  kstri,
							  fid);

		else if (!strncmp(kstri,"<",1))
			KML_Object::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for (ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
