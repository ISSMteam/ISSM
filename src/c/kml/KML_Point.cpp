/*!\file KML_Point.cpp
 * \brief: implementation of the kml_point object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Point.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Point::KML_Point(){/*{{{*/

	extrude   =false;
	memcpy(altmode,"clampToGround",(strlen("clampToGround")+1)*sizeof(char));

	coords[0] = 0.;
	coords[1] = 0.;
	coords[2] = 0.;

}
/*}}}*/
KML_Point::~KML_Point(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_Point::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Point:\n");
	KML_Geometry::Echo();

	if(flag) _printf0_("       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_("       altmode: \"" << altmode << "\"\n");
	if(flag) _printf0_("        coords: (" << coords[0] << "," << coords[1] << "," << coords[2] << ")\n");

	return;
}
/*}}}*/
void  KML_Point::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Point::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Point::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Point:\n");
	KML_Geometry::DeepEcho(indent);

	if(flag) _printf0_(indent << "       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "       altmode: \"" << altmode << "\"\n");
	if(flag) _printf0_(indent << "        coords: (" << coords[0] << "," << coords[1] << "," << coords[2] << ")\n");

	return;
}
/*}}}*/
void  KML_Point::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<Point",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Geometry::Write(filout,indent);

	fprintf(filout,"%s  <extrude>%d</extrude>\n",indent,(extrude ? 1 : 0));
	fprintf(filout,"%s  <altitudeMode>%s</altitudeMode>\n",indent,altmode);
	fprintf(filout,"%s  <coordinates>%0.16g,%0.16g,%0.16g</coordinates>\n",
			indent,coords[0],coords[1],coords[2]);

	fprintf(filout,"%s</Point>\n",indent);

	return;
}
/*}}}*/
void  KML_Point::Read(FILE* fid,char* kstr){/*{{{*/

	double*      pcoords=&coords[0];
	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Point", 7)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Point::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Point::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<extrude>"))
			KMLFileTokenParse(&extrude   , kstri, fid);
		else if (!strcmp(kstri,"<altitudeMode>"))
			KMLFileTokenParse( altmode   ,NULL,KML_POINT_ALTMODE_LENGTH, kstri, fid);
		else if (!strcmp(kstri,"<coordinates>"))
			KMLFileTokenParse(&pcoords   ,NULL,3, kstri, fid);

		else if (!strncmp(kstri,"<",1))
			KML_Geometry::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
void  KML_Point::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int     i;
	double  lat,lon,x,y;
	char    nstr2[81];

/*  extract latitude and longitude  */

	lon=coords[0];
	lat=coords[1];

/*  convert latitude and longitude to x and y  */

	if (sgn) {
		Ll2xyx(&x,&y,&lat,&lon,1,sgn,cm,sp);
	}
	else {
		memcpy(&x,&lon,1*sizeof(IssmDouble));
		memcpy(&y,&lat,1*sizeof(IssmDouble));
	}

/*  write header  */

	memcpy(nstr2,nstr,(strlen(nstr)+1)*sizeof(char));

	for (i=0; i<strlen(nstr2); i++)
		if ((nstr2[i] == ' ') || (nstr2[i] == '\t'))
			nstr2[i]='_';
	fprintf(fid,"## Name:%s\n",nstr2);
	fprintf(fid,"## Icon:0\n");
	fprintf(fid,"# Points Count	Value\n");
    fprintf(fid,"%u	%s\n",1,"1.");
	fprintf(fid,"# X pos	Y pos\n");

/*  write vertex  */

    fprintf(fid,"%lf\t%lf\n",x,y);

/*  write blank line  */

	fprintf(fid,"\n");

	return;
}
/*}}}*/
