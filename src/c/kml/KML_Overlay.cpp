/*!\file KML_Overlay.cpp
 * \brief: implementation of the kml_overlay abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KMLFileReadUtils.h"
#include "./KML_Overlay.h"
#include "./KML_Icon.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Overlay::KML_Overlay(){/*{{{*/

	strcpy(color     ,"ffffffff");
	memcpy(color,"ffffffff",(strlen("ffffffff")+1)*sizeof(char));

	draword   = 0;
	icon      =NULL;

}
/*}}}*/
KML_Overlay::~KML_Overlay(){/*{{{*/

	if (icon) {
		delete icon;
		icon      =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Overlay::Echo(){/*{{{*/

	KML_Feature::Echo();
	_printf0_("         color: \"" << color << "\"\n");
	_printf0_("       draword: " << draword << "\n");
	_printf0_("          icon: " << icon << "\n");
}
/*}}}*/
void  KML_Overlay::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Overlay::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Overlay::DeepEcho(const char* indent){/*{{{*/

	char  indent2[81];
	KML_Feature::DeepEcho(indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	_printf0_(indent << "         color: " << color << "\n");
	_printf0_(indent << "       draword: " << draword << "\n");
	if (icon)
		icon->DeepEcho(indent2);
	else
		_printf0_(indent << "          icon: " << icon << "\n");
}
/*}}}*/
void  KML_Overlay::Write(FILE* filout,const char* indent){/*{{{*/

	char  indent2[81];

	KML_Feature::Write(filout,indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

	if (color     && strlen(color))
		fprintf(filout,"%s  <color>%s</color>\n",indent,color);
	fprintf(filout,"%s  <drawOrder>%d</drawOrder>\n",indent,draword);
	if (icon)
		icon->Write(filout,indent2);

	return;
}
/*}}}*/
void  KML_Overlay::Read(FILE* fid,char* kstr){/*{{{*/

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</Overlay", 9)) {
		xDelete<char>(kstr);
		return;
	}
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_Overlay::Read -- Unexpected closing tag " << kstr << ".\n");}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_Overlay::Read -- Unexpected field \"" << kstr << "\".\n");}

	else if (!strcmp(kstr,"<color>"))
		KMLFileTokenParse( color     ,NULL,KML_OVERLAY_COLOR_LENGTH,
						  kstr,
						  fid);
	else if (!strcmp(kstr,"<drawOrder>"))
		KMLFileTokenParse(&draword   ,
						  kstr,
						  fid);

	else if (!strncmp(kstr,"<Icon", 5)) {
		icon      =new KML_Icon();
		icon      ->Read(fid,kstr);
	}

	else if (!strncmp(kstr,"<",1))
		KML_Feature::Read(fid,kstr);

	return;
}
/*}}}*/
