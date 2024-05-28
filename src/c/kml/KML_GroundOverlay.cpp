/*!\file KML_GroundOverlay.cpp
 * \brief: implementation of the kml_groundoverlay object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Object.h"
#include "./KML_LatLonBox.h"
#include "./KML_GroundOverlay.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_GroundOverlay::KML_GroundOverlay(){/*{{{*/

	altitude  = 0.;
	memcpy(altmode,"clampToGround",(strlen("clampToGround")+1)*sizeof(char));

	llbox     =NULL;

}
/*}}}*/
KML_GroundOverlay::~KML_GroundOverlay(){/*{{{*/

	if (llbox) {
		delete llbox;
		llbox     =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_GroundOverlay::Echo(){/*{{{*/

	_printf_("KML_GroundOverlay:\n");
	KML_Overlay::Echo();

	_printf_("         altitude: " << altitude << "\n");
	_printf_("          altmode: " << altmode << "\n");
	_printf_("            llbox: " << llbox << "\n");
}
/*}}}*/
void  KML_GroundOverlay::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_GroundOverlay::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_GroundOverlay::DeepEcho(const char* indent){/*{{{*/

	char  indent2[81];

	_printf_(indent << "KML_GroundOverlay:\n");
	KML_Overlay::DeepEcho(indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	_printf_(indent<<"      altitude: " << altitude << "\n");
	_printf_(indent<<"       altmode: " << altmode << "\n");
	if (llbox)
	 llbox->DeepEcho(indent2);
	else
	 _printf_(indent<<"         llbox: " << llbox << "\n");
}
/*}}}*/
void  KML_GroundOverlay::Write(FILE* filout,const char* indent){/*{{{*/

	char  indent2[81];

	fprintf(filout,"%s<GroundOverlay",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Overlay::Write(filout,indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

	fprintf(filout,"%s  <altitude>%0.16g</altitude>\n",indent,altitude);
	fprintf(filout,"%s  <altitudeMode>%s</altitudeMode>\n",indent,altmode);
	if (llbox)
		llbox->Write(filout,indent2);

	fprintf(filout,"%s</GroundOverlay>\n",indent);

	return;
}
/*}}}*/
void  KML_GroundOverlay::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</GroundOverlay",15)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_GroundOverlay::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_GroundOverlay::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<altitude>"))
			KMLFileTokenParse(&altitude  ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<altitudeMode>"))
			KMLFileTokenParse( altmode   ,NULL,KML_GROUNDOVERLAY_ALTMODE_LENGTH,
							  kstri,
							  fid);
		else if (!strncmp(kstri,"<LatLonBox",10)) {
			llbox     =new KML_LatLonBox();
			llbox     ->Read(fid,kstri);
		}

		else if (!strncmp(kstri,"<",1))
			KML_Overlay::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
