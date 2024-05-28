/*!\file KML_Icon.cpp
 * \brief: implementation of the kml_feature abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Icon.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Icon::KML_Icon(){/*{{{*/

	strcpy(href      ,"");
	strcpy(refmode   ,"onChange");
	refint    = 4.;
	strcpy(vrefmode  ,"never");
	vreftime  = 4.;
	vboundsc  = 1.;
	strcpy(vformat   ,"");
	strcpy(hquery    ,"");

}
/*}}}*/
KML_Icon::~KML_Icon(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_Icon::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Icon:\n");
	KML_Object::Echo();

	if(flag) _printf0_("          href: \"" << href << "\"\n");
	if(flag) _printf0_("       refmode: \"" << refmode << "\"\n");
	if(flag) _printf0_("        refint: " << refint << "\n");
	if(flag) _printf0_("      vrefmode: \"" << vrefmode << "\"\n");
	if(flag) _printf0_("      vreftime: " << vreftime << "\n");
	if(flag) _printf0_("      vboundsc: " << vboundsc << "\n");
	if(flag) _printf0_("       vformat: \"" << vformat << "\"\n");
	if(flag) _printf0_("        hquery: \"" << hquery << "\"\n");

	return;
}
/*}}}*/
void  KML_Icon::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Icon::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Icon::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Icon:\n");
	KML_Object::DeepEcho(indent);

	if(flag) _printf0_(indent << "          href: \"" << href << "\"\n");
	if(flag) _printf0_(indent << "       refmode: \"" << refmode << "\"\n");
	if(flag) _printf0_(indent << "        refint: " << refint << "\n");
	if(flag) _printf0_(indent << "      vrefmode: \"" << vrefmode << "\"\n");
	if(flag) _printf0_(indent << "      vreftime: " << vreftime << "\n");
	if(flag) _printf0_(indent << "      vboundsc: " << vboundsc << "\n");
	if(flag) _printf0_(indent << "       vformat: \"" << vformat << "\"\n");
	if(flag) _printf0_(indent << "        hquery: \"" << hquery << "\"\n");

	return;
}
/*}}}*/
void  KML_Icon::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<Icon",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Object::Write(filout,indent);

	if (href     && strlen(href))
		fprintf(filout,"%s  <href>%s</href>\n",indent,href);
	if (refmode  && strlen(refmode))
		fprintf(filout,"%s  <refreshMode>%s</refreshMode>\n",indent,refmode);
	fprintf(filout,"%s  <refreshInterval>%g</refreshInterval>\n",indent,refint);
	if (vrefmode && strlen(vrefmode))
		fprintf(filout,"%s  <viewRefreshMode>%s</viewRefreshMode>\n",indent,vrefmode);
	fprintf(filout,"%s  <viewRefreshTime>%g</viewRefreshTime>\n",indent,vreftime);
	fprintf(filout,"%s  <viewBoundScale>%g</viewBoundScale>\n",indent,vboundsc);
	if (vformat  && strlen(vformat))
		fprintf(filout,"%s  <viewFormat>%s</viewFormat>\n",indent,vformat);
	if (hquery   && strlen(hquery))
		fprintf(filout,"%s  <httpQuery>%s</httpQuery>\n",indent,hquery);

	fprintf(filout,"%s</Icon>\n",indent);

	return;
}
/*}}}*/
void  KML_Icon::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Icon", 6)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Icon::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Icon::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<href>"))
			KMLFileTokenParse( href      ,NULL,KML_ICON_HREF_LENGTH, kstri, fid);
		else if (!strcmp(kstri,"<refreshMode>"))
			KMLFileTokenParse( refmode   ,NULL,KML_ICON_REFMODE_LENGTH, kstri, fid);
		else if (!strcmp(kstri,"<refreshInterval>"))
			KMLFileTokenParse(&refint    , kstri, fid);
		else if (!strcmp(kstri,"<viewRefreshMode>"))
			KMLFileTokenParse( vrefmode  ,NULL,KML_ICON_VREFMODE_LENGTH, kstri, fid);
		else if (!strcmp(kstri,"<viewRefreshTime>"))
			KMLFileTokenParse(&vreftime  , kstri, fid);
		else if (!strcmp(kstri,"<viewBoundScale>"))
			KMLFileTokenParse(&vboundsc  , kstri, fid);
		else if (!strcmp(kstri,"<viewFormat>"))
			KMLFileTokenParse( vformat   ,NULL,KML_ICON_VFORMAT_LENGTH, kstri, fid);
		else if (!strcmp(kstri,"<httpQuery>"))
			KMLFileTokenParse( hquery    ,NULL,KML_ICON_HQUERY_LENGTH, kstri, fid);

		else if (!strncmp(kstri,"<",1))
			KML_Object::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
