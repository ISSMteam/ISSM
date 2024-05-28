/*!\file KML_MultiGeometry.cpp
 * \brief: implementation of the kml_multigeometry object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Object.h"
#include "./KML_Point.h"
#include "./KML_Polygon.h"
#include "./KML_LineString.h"
#include "./KML_LinearRing.h"
#include "./KMLFileReadUtils.h"
#include "./KML_MultiGeometry.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_MultiGeometry::KML_MultiGeometry(){/*{{{*/

	geometry  =new DataSet;

}
/*}}}*/
KML_MultiGeometry::~KML_MultiGeometry(){/*{{{*/

	if (geometry) {
		delete geometry;
		geometry  =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_MultiGeometry::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Multigeometry:\n");
	KML_Geometry::Echo();

	if(flag) _printf0_("      geometry: (size=" << geometry->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_MultiGeometry::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_MultiGeometry::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_MultiGeometry::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Multigeometry:\n");
	KML_Geometry::DeepEcho(indent);

/*  loop over the geometry elements for the multigeometry  */

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
void  KML_MultiGeometry::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent2[81];

	fprintf(filout,"%s<MultiGeometry",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Geometry::Write(filout,indent);

/*  loop over the geometry elements for the multigeometry  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

	for (i=0; i<geometry->Size(); i++)
		((KML_Geometry *)geometry->GetObjectByOffset(i))->Write(filout,indent2);

	fprintf(filout,"%s</MultiGeometry>\n",indent);

	return;
}
/*}}}*/
void  KML_MultiGeometry::Read(FILE* fid,char* kstr){/*{{{*/

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
		if      (!strncmp(kstri,"</MultiGeometry",15)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_MultiGeometry::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_MultiGeometry::Read -- Unexpected field \"" << kstri << "\".\n");}

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
void  KML_MultiGeometry::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int   i;

/*  loop over the geometry elements for the multigeometry  */

	for (i=0; i<geometry->Size(); i++)
		((KML_Object *)geometry->GetObjectByOffset(i))->WriteExp(fid,nstr,sgn,cm,sp);

	return;
}
/*}}}*/
