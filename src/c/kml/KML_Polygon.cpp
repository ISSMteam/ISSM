/*!\file KML_Polygon.cpp
 * \brief: implementation of the kml_polygon object
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
#include "./KML_Polygon.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Polygon::KML_Polygon(){/*{{{*/

	extrude   =false;
	tessellate=false;
	memcpy(altmode,"clampToGround",(strlen("clampToGround")+1)*sizeof(char));

	outer     =new DataSet;
	inner     =new DataSet;

}
/*}}}*/
KML_Polygon::~KML_Polygon(){/*{{{*/

	if (inner) {
		delete inner;
		inner     =NULL;
	}

	if (outer) {
		delete outer;
		outer     =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Polygon::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Polygon:\n");
	KML_Geometry::Echo();

	if(flag) _printf0_("       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_("    tessellate: " << (tessellate ? "true" : "false") << "\n");
	if(flag) _printf0_("       altmode: \"" << altmode << "\"\n");
	if(flag) _printf0_("         outer: (size=" << outer->Size() << ")\n");
	if(flag) _printf0_("         inner: (size=" << inner->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_Polygon::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Polygon::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Polygon::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Polygon:\n");
	KML_Geometry::DeepEcho(indent);

	if(flag) _printf0_(indent << "       extrude: " << (extrude ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "    tessellate: " << (tessellate ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "       altmode: \"" << altmode << "\"\n");

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (outer->Size())
		for (i=0; i<outer->Size(); i++) {
			if(flag) _printf0_(indent << "         outer: -------- begin [" << i << "] --------\n");
			((KML_LinearRing *)outer->GetObjectByOffset(i))->DeepEcho(indent2);
			if(flag) _printf0_(indent << "         outer: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "         outer: [empty]\n");

	if (inner->Size())
		for (i=0; i<inner->Size(); i++) {
			if(flag) _printf0_(indent << "         inner: -------- begin [" << i << "] --------\n");
			((KML_LinearRing *)inner->GetObjectByOffset(i))->DeepEcho(indent2);
			if(flag) _printf0_(indent << "         inner: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "         inner: [empty]\n");

	return;
}
/*}}}*/
void  KML_Polygon::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent4[81];

	fprintf(filout,"%s<Polygon",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Geometry::Write(filout,indent);

	fprintf(filout,"%s  <extrude>%d</extrude>\n",indent,(extrude ? 1 : 0));
	fprintf(filout,"%s  <tessellate>%d</tessellate>\n",indent,(tessellate ? 1 : 0));
	fprintf(filout,"%s  <altitudeMode>%s</altitudeMode>\n",indent,altmode);

	memcpy(indent4,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent4,"    ");

/*  check outer boundary for the polygon  */

	fprintf(filout,"%s  <outerBoundaryIs>\n",indent);
	if (outer->Size())
		((KML_LinearRing *)outer->GetObjectByOffset(0))->Write(filout,indent4);
	fprintf(filout,"%s  </outerBoundaryIs>\n",indent);

/*  loop over any inner boundaries for the polygon  */

	for (i=0; i<inner->Size(); i++) {
		fprintf(filout,"%s  <innerBoundaryIs>\n",indent);
		((KML_LinearRing *)inner->GetObjectByOffset(i))->Write(filout,indent4);
		fprintf(filout,"%s  </innerBoundaryIs>\n",indent);
	}

	fprintf(filout,"%s</Polygon>\n",indent);

	return;
}
/*}}}*/
void  KML_Polygon::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	char*        kstrj;
	int          ncom=0;
	char**       pcom=NULL;
	KML_Object*  kobj;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Polygon", 9)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Polygon::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Polygon::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strcmp(kstri,"<extrude>"))
			KMLFileTokenParse(&extrude   ,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<tessellate>"))
			KMLFileTokenParse(&tessellate,
							  kstri,
							  fid);
		else if (!strcmp(kstri,"<altitudeMode>"))
			KMLFileTokenParse( altmode   ,NULL,KML_POLYGON_ALTMODE_LENGTH,
							  kstri,
							  fid);

		else if (!strcmp(kstri,"<outerBoundaryIs>"))

/*  loop over and process fields within outer boundary  */

			while((kstrj=KMLFileToken(fid, &ncom,&pcom))){
				if      (!strncmp(kstrj,"</outerBoundaryIs",17)) {
					xDelete<char>(kstrj);
					break;
				}
				else if (!strncmp(kstrj,"</",2))
				  {_error_("KML_Polygon::Read -- Unexpected closing tag " << kstrj << ".\n");}
				else if (strncmp(kstrj,"<",1))
				  {_error_("KML_Polygon::Read -- Unexpected field \"" << kstrj << "\".\n");}

				else if (!strncmp(kstrj,"<LinearRing",11)) {
					kobj=(KML_Object*)new KML_LinearRing();
					kobj->Read(fid,kstrj);
					outer     ->AddObject((Object*)kobj);
				}

				else if (!strncmp(kstrj,"<",1))
					KML_Geometry::Read(fid,kstrj);

				xDelete<char>(kstrj);
			}

		else if (!strcmp(kstri,"<innerBoundaryIs>"))

/*  loop over and process fields within inner boundaries  */

			while((kstrj=KMLFileToken(fid, &ncom,&pcom))){
				if      (!strncmp(kstrj,"</innerBoundaryIs",17)) {
					xDelete<char>(kstrj);
					break;
				}
				else if (!strncmp(kstrj,"</",2))
				  {_error_("KML_Polygon::Read -- Unexpected closing tag " << kstrj << ".\n");}
				else if (strncmp(kstrj,"<",1))
				  {_error_("KML_Polygon::Read -- Unexpected field \"" << kstrj << "\".\n");}

				else if (!strncmp(kstrj,"<LinearRing",11)) {
					kobj=(KML_Object*)new KML_LinearRing();
					kobj->Read(fid,kstrj);
					inner     ->AddObject((Object*)kobj);
				}

				else if (!strncmp(kstrj,"<",1))
					KML_Geometry::Read(fid,kstrj);

				xDelete<char>(kstrj);
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
void  KML_Polygon::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int   i;
	char  nstr2[81];

/*  check outer boundary for the polygon  */

	if (outer->Size()) {
		if (strlen(nstr))
			sprintf(nstr2,"%s (outer)",nstr);
		else
			sprintf(nstr2,"(outer)");

		((KML_LinearRing *)outer->GetObjectByOffset(0))->WriteExp(fid,nstr2,sgn,cm,sp);
	}

/*  loop over any inner boundaries for the polygon  */

	for (i=0; i<inner->Size(); i++) {
		if (strlen(nstr))
			sprintf(nstr2,"%s (inner %d of %d)",nstr,i+1,inner->Size());
		else
			sprintf(nstr2,"(inner %d of %d)",i+1,inner->Size());

		((KML_LinearRing *)inner->GetObjectByOffset(i))->WriteExp(fid,nstr2,sgn,cm,sp);
	}

	return;
}
/*}}}*/
