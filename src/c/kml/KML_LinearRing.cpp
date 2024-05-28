/*!\file KML_LinearRing.cpp
 * \brief: implementation of the kml_linearring object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_LinearRing.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_LinearRing::KML_LinearRing(){/*{{{*/

	extrude   =false;
	tessellate=false;
	memcpy(altmode,"clampToGround",(strlen("clampToGround")+1)*sizeof(char));

	ncoord    =0;
	coords    =NULL;

}
/*}}}*/
KML_LinearRing::~KML_LinearRing(){/*{{{*/

	if (coords) xDelete<double>(coords);

	coords    =NULL;
	ncoord    =0;

}
/*}}}*/

/*Other*/
void  KML_LinearRing::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_LinearRing:\n");
	KML_Geometry::Echo();

	if(flag) _printf0_("       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_("    tessellate: " << (tessellate ? "true" : "false") << "\n");
	if(flag) _printf0_("       altmode: \"" << altmode << "\"\n");
	if(flag) _printf0_("        coords: (ncoord=" << ncoord << ")\n");

	return;
}
/*}}}*/
void  KML_LinearRing::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_LinearRing::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_LinearRing::DeepEcho(const char* indent){/*{{{*/

	int   i;
	bool  flag=true;

	if(flag) _printf0_(indent << "KML_LinearRing:\n");
	KML_Geometry::DeepEcho(indent);

	if(flag) _printf0_(indent << "       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "    tessellate: " << (tessellate ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "       altmode: \"" << altmode << "\"\n");
	if(flag) _printf0_(indent << "        coords: (ncoord=" << ncoord << ")\n");
	for (i=0; i<ncoord; i++)
		if(flag)_printf_(indent << "                (" <<coords[3*i+0] << "," <<coords[3*i+1] << "," <<coords[3*i+2] << ")\n\n");

	return;
}
/*}}}*/
void  KML_LinearRing::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;

	fprintf(filout,"%s<LinearRing",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Geometry::Write(filout,indent);

	fprintf(filout,"%s  <extrude>%d</extrude>\n",indent,(extrude ? 1 : 0));
	fprintf(filout,"%s  <tessellate>%d</tessellate>\n",indent,(tessellate ? 1 : 0));
	fprintf(filout,"%s  <altitudeMode>%s</altitudeMode>\n",indent,altmode);
	fprintf(filout,"%s  <coordinates>\n",indent);

/*  loop over the coordinates for the linearring  */

	for (i=0; i<ncoord; i++)
		fprintf(filout,"%s    %0.16g,%0.16g,%0.16g\n",indent,coords[3*i+0],coords[3*i+1],coords[3*i+2]);

	fprintf(filout,"%s  </coordinates>\n",indent);
	fprintf(filout,"%s</LinearRing>\n",indent);

	return;
}
/*}}}*/
void  KML_LinearRing::Read(FILE* fid,char* kstr){/*{{{*/

	char  *kstri = NULL;
	int    ncom  = 0;
	char **pcom  = NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,kstr)) return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid,&ncom,&pcom))){
		if (!strncmp(kstri,"</LinearRing",12)){
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_LinearRing::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_LinearRing::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<extrude>"))
			KMLFileTokenParse(&extrude,kstri,fid);
		else if (!strcmp(kstri,"<tessellate>"))
			KMLFileTokenParse(&tessellate,kstri,fid);
		else if (!strcmp(kstri,"<altitudeMode>"))
			KMLFileTokenParse(altmode,NULL,KML_LINEARRING_ALTMODE_LENGTH,kstri,fid);
		else if (!strcmp(kstri,"<coordinates>"))
			KMLFileTokenParse(&coords,&ncoord,3,0,kstri,fid);
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
void  KML_LinearRing::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int     i;
	double  *lat,*lon,*x,*y;
	char    nstr2[81];

/*  extract latitude and longitude into vectors  */

	lat=xNew<IssmPDouble>(ncoord);
	lon=xNew<IssmPDouble>(ncoord);
	for (i=0; i<ncoord; i++) {
		lon[i]=coords[3*i+0];
		lat[i]=coords[3*i+1];
	}

/*  convert latitude and longitude to x and y  */

	x  =xNew<IssmPDouble>(ncoord);
	y  =xNew<IssmPDouble>(ncoord);
	if (sgn) {
		Ll2xyx(x,y,lat,lon,ncoord,sgn,cm,sp);
	}
	else {
		memcpy(x,lon,ncoord*sizeof(IssmDouble));
		memcpy(y,lat,ncoord*sizeof(IssmDouble));
	}

/*  write header  */

	memcpy(nstr2,nstr,(strlen(nstr)+1)*sizeof(char));

	for (i=0; i<strlen(nstr2); i++)
		if ((nstr2[i] == ' ') || (nstr2[i] == '\t'))
			nstr2[i]='_';
	fprintf(fid,"## Name:%s\n",nstr2);
	fprintf(fid,"## Icon:0\n");
	fprintf(fid,"# Points Count	Value\n");
	if ((lat[ncoord-1] != lat[0]) || (lon[ncoord-1] != lon[0]))
	    fprintf(fid,"%u	%s\n",ncoord+1,"1.");
	else
	    fprintf(fid,"%u	%s\n",ncoord  ,"1.");
	fprintf(fid,"# X pos	Y pos\n");

/*  write vertices, making sure ring is closed  */

	for (i=0; i<ncoord; i++)
	    fprintf(fid,"%lf\t%lf\n",x[i],y[i]);
	if ((lat[ncoord-1] != lat[0]) || (lon[ncoord-1] != lon[0]))
	    fprintf(fid,"%lf\t%lf\n",x[0],y[0]);

/*  write blank line  */

	fprintf(fid,"\n");

	xDelete<double>(y);
	xDelete<double>(x);
	xDelete<double>(lon);
	xDelete<double>(lat);

	return;
}
/*}}}*/
