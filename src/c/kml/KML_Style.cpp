/*!\file KML_Style.cpp
 * \brief: implementation of the kml_style object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KMLFileReadUtils.h"
#include "./KML_LineStyle.h"
#include "./KML_PolyStyle.h"
#include "./KML_Style.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Style::KML_Style(){/*{{{*/

	icon      =NULL;
	label     =NULL;
	line      =NULL;
	poly      =NULL;
	balloon   =NULL;
	list      =NULL;

}
/*}}}*/
KML_Style::~KML_Style(){/*{{{*/

	if (list) {
//		delete list;
		list      =NULL;
	}
	if (balloon) {
//		delete balloon;
		balloon   =NULL;
	}
	if (poly) {
		delete poly;
		poly      =NULL;
	}
	if (line) {
		delete line;
		line      =NULL;
	}
	if (label) {
//		delete label;
		label     =NULL;
	}
	if (icon) {
//		delete icon;
		icon      =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Style::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Style:\n");
	KML_StyleSelector::Echo();

	if(flag) _printf0_("          icon: " << icon << "\n");
	if(flag) _printf0_("         label: " << label << "\n");
	if(flag) _printf0_("          line: " << line << "\n");
	if(flag) _printf0_("          poly: " << poly << "\n");
	if(flag) _printf0_("       balloon: " << balloon << "\n");
	if(flag) _printf0_("          list: " << list << "\n");

	return;
}
/*}}}*/
void  KML_Style::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Style::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Style::DeepEcho(const char* indent){/*{{{*/

	char  indent2[81];
	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Style:\n");
	KML_StyleSelector::DeepEcho(indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

//	if (icon)
//		icon->DeepEcho(indent2);
//	else
		if(flag) _printf0_(indent << "          icon: " << icon << "\n");
//	if (label)
//		label->DeepEcho(indent2);
//	else
		if(flag) _printf0_(indent << "         label: " << label << "\n");
	if (line)
		line->DeepEcho(indent2);
	else
		if(flag) _printf0_(indent << "          line: " << line << "\n");
	if (poly)
		poly->DeepEcho(indent2);
	else
		if(flag) _printf0_(indent << "          poly: " << poly << "\n");
//	if (balloon)
//		balloon->DeepEcho(indent2);
//	else
		if(flag) _printf0_(indent << "       balloon: " << balloon << "\n");
//	if (list)
//		list->DeepEcho(indent2);
//	else
		if(flag) _printf0_(indent << "          list: " << list << "\n");

	return;
}
/*}}}*/
void  KML_Style::Write(FILE* filout,const char* indent){/*{{{*/

	char  indent2[81];

	fprintf(filout,"%s<Style",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_StyleSelector::Write(filout,indent);

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

//	if (icon)
//		icon->Write(filout,indent2);
//	if (label)
//		label->Write(filout,indent2);
	if (line)
		line->Write(filout,indent2);
	if (poly)
		poly->Write(filout,indent2);
//	if (balloon)
//		balloon->Write(filout,indent2);
//	if (list)
//		list->Write(filout,indent2);

	fprintf(filout,"%s</Style>\n",indent);

	return;
}
/*}}}*/
void  KML_Style::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Style", 7)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Style::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Style::Read -- Unexpected field \"" << kstri << "\".\n");}

//		else if (!strncmp(kstri,"<IconStyle",10)) {
//			icon      =new KML_IconStyle();
//			icon      ->Read(fid,kstri);
//		}

//		else if (!strncmp(kstri,"<LabelStyle",11)) {
//			label     =new KML_LabelStyle();
//			label     ->Read(fid,kstri);
//		}

		else if (!strncmp(kstri,"<LineStyle",10)) {
			line      =new KML_LineStyle();
			line      ->Read(fid,kstri);
		}

		else if (!strncmp(kstri,"<PolyStyle",10)) {
			poly      =new KML_PolyStyle();
			poly      ->Read(fid,kstri);
		}

//		else if (!strncmp(kstri,"<BalloonStyle",13)) {
//			balloon   =new KML_BalloonStyle();
//			balloon   ->Read(fid,kstri);
//		}

//		else if (!strncmp(kstri,"<ListStyle",10)) {
//			list      =new KML_ListStyle();
//			list      ->Read(fid,kstri);
//		}

		else if (!strncmp(kstri,"<",1))
			KML_StyleSelector::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
