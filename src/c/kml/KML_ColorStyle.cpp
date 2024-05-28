/*!\file KML_ColorStyle.cpp
 * \brief: implementation of the kml_colorstyle abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_ColorStyle.h"
#include "./KML_SubStyle.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_ColorStyle::KML_ColorStyle(){/*{{{*/

	strcpy(color     ,"ffffffff");
	strcpy(colormode ,"normal");

}
/*}}}*/
KML_ColorStyle::~KML_ColorStyle(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_ColorStyle::Echo(){/*{{{*/

	bool  flag=true;

	KML_SubStyle::Echo();

	if(flag) _printf0_("         color: " << color << "\n");
	if(flag) _printf0_("     colormode: " << colormode << "\n");

	return;
}
/*}}}*/
void  KML_ColorStyle::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_ColorStyle::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_ColorStyle::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	KML_SubStyle::DeepEcho(indent);

	if(flag) _printf0_(indent << "         color: " << color << "\n");
	if(flag) _printf0_(indent << "     colormode: " << colormode << "\n");
}
/*}}}*/
void  KML_ColorStyle::Write(FILE* filout,const char* indent){/*{{{*/

	KML_SubStyle::Write(filout,indent);

	if (color     && strlen(color))
		fprintf(filout,"%s  <color>%s</color>\n",indent,color);
	if (colormode && strlen(colormode))
		fprintf(filout,"%s  <colorMode>%s</colorMode>\n",indent,colormode);

	return;
}
/*}}}*/
void  KML_ColorStyle::Read(FILE* fid,char* kstr){/*{{{*/

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</ColorStyle",12))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_ColorStyle::Read -- Unexpected closing tag " << kstr);}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_ColorStyle::Read -- Unexpected field \"" << kstr << "\"");}

	else if (!strcmp(kstr,"<color>"))
		KMLFileTokenParse( color     ,NULL,KML_COLORSTYLE_COLOR_LENGTH, kstr, fid);
	else if (!strcmp(kstr,"<colorMode>"))
		KMLFileTokenParse( colormode ,NULL,KML_COLORSTYLE_COLORMODE_LENGTH, kstr, fid);

	else if (!strncmp(kstr,"<",1))
		KML_SubStyle::Read(fid,kstr);

	return;
}
/*}}}*/
