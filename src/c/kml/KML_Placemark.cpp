/*!\file KML_Placemark.cpp
 * \brief: implementation of the kml_placemark object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KMLFileReadUtils.h"
#include "./KML_Geometry.h"
#include "./KML_Point.h"
#include "./KML_LineString.h"
#include "./KML_Polygon.h"
#include "./KML_MultiGeometry.h"
#include "./KML_LinearRing.h"
#include "./KML_Placemark.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Placemark::KML_Placemark(){/*{{{*/

	geometry  =new DataSet;

}
/*}}}*/
KML_Placemark::~KML_Placemark(){/*{{{*/

	if (geometry) {
		delete geometry;
		geometry  =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Placemark::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Placemark:\n");
	KML_Feature::Echo();

	if(flag) _printf0_("      geometry: (size=" << geometry->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_Placemark::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Placemark::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Placemark::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Placemark:\n");
	KML_Feature::DeepEcho(indent);

/*  loop over the geometry elements for the placemark  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (geometry->Size())
		for (i=0; i<geometry->Size(); i++) {
			if(flag) _printf0_(indent << "      geometry: -------- begin [" << i << "] --------\n");
			((KML_Geometry *)geometry->GetObjectByOffset(i))->DeepEcho(indent2);
			if(flag) _printf0_(indent << "      geometry: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "      geometry: [empty]\n");

	return;
}
/*}}}*/
void  KML_Placemark::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent2[81];

	fprintf(filout,"%s<Placemark",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Feature::Write(filout,indent);

/*  loop over the geometry elements for the placemark  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

	for (i=0; i<geometry->Size(); i++)
		((KML_Geometry *)geometry->GetObjectByOffset(i))->Write(filout,indent2);

	fprintf(filout,"%s</Placemark>\n",indent);

	return;
}
/*}}}*/
void  KML_Placemark::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;
	KML_Object*  kobj;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Placemark",11)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Placemark::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Placemark::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strncmp(kstri,"<Point", 6)) {
			kobj=(KML_Object*)new KML_Point();
			kobj->Read(fid,kstri);
			geometry  ->AddObject((Object*)kobj);
		}

		else if (!strncmp(kstri,"<LineString",11)) {
			kobj=(KML_Object*)new KML_LineString();
			kobj->Read(fid,kstri);
			geometry  ->AddObject((Object*)kobj);
		}

		else if (!strncmp(kstri,"<LinearRing",11)) {
			kobj=(KML_Object*)new KML_LinearRing();
			kobj->Read(fid,kstri);
			geometry  ->AddObject((Object*)kobj);
		}

		else if (!strncmp(kstri,"<Polygon", 8)) {
			kobj=(KML_Object*)new KML_Polygon();
			kobj->Read(fid,kstri);
			geometry  ->AddObject((Object*)kobj);
		}

		else if (!strncmp(kstri,"<MultiGeometry",14)) {
			kobj=(KML_Object*)new KML_MultiGeometry();
			kobj->Read(fid,kstri);
			geometry  ->AddObject((Object*)kobj);
		}

		else if (!strncmp(kstri,"<",1))
			KML_Feature::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
void  KML_Placemark::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int   i;
	char  nstr2[81];

/*  loop over the geometry elements for the placemark  */

	for (i=0; i<geometry->Size(); i++) {
		if (strlen(nstr))
			sprintf(nstr2,"%s %s",nstr,name);
		else
			sprintf(nstr2,"%s",name);

		((KML_Object *)geometry->GetObjectByOffset(i))->WriteExp(fid,nstr2,sgn,cm,sp);
	}

	return;
}
/*}}}*/
