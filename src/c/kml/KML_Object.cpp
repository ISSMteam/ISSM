/*!\file KML_Object.cpp
 * \brief: implementation of the kml_object abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Object.h"
#include "./KML_Attribute.h"
#include "./KML_Comment.h"
#include "./KML_Unknown.h"
#include "./KML_LatLonBox.h"
#include "./KML_Icon.h"
#include "./KML_MultiGeometry.h"
#include "./KML_Document.h"
#include "./KML_LinearRing.h"
#include "./KML_LineStyle.h"
#include "./KML_LineString.h"
#include "./KML_PolyStyle.h"
#include "./KML_Polygon.h"
#include "./KML_Point.h"
#include "./KML_GroundOverlay.h"
#include "./KML_Placemark.h"
#include "./KML_Folder.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Object::KML_Object(){/*{{{*/

	attrib    =new DataSet;
	commnt    =new DataSet;
	kmlobj    =new DataSet;

}
/*}}}*/
KML_Object::~KML_Object(){/*{{{*/

	if (attrib) {
		delete attrib;
		attrib    =NULL;
	}
	if (commnt) {
		delete commnt;
		commnt    =NULL;
	}
	if (kmlobj) {
		delete kmlobj;
		kmlobj    =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Object::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("        attrib: (size=" << attrib->Size() << ")\n");
	if(flag) _printf0_("        commnt: (size=" << commnt->Size() << ")\n");
	if(flag) _printf0_("        kmlobj: (size=" << kmlobj->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_Object::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Object::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Object::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

/*  loop over the attributes for the object  */

	if (attrib->Size())
		for (i=0; i<attrib->Size(); i++) {
			((KML_Attribute *)attrib->GetObjectByOffset(i))->DeepEcho(indent);
		}
	else
		if(flag) _printf0_(indent << "        attrib: [empty]\n");

/*  loop over the comments for the object  */

	if (commnt->Size())
		for (i=0; i<commnt->Size(); i++) {
			((KML_Comment *)commnt->GetObjectByOffset(i))->DeepEcho(indent);
		}
	else
		if(flag) _printf0_(indent << "        commnt: [empty]\n");

/*  loop over the unknown objects for the object  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (kmlobj->Size())
		for (i=0; i<kmlobj->Size(); i++) {
            if(flag) _printf0_(indent << "        kmlobj: -------- begin [" << i << "] --------\n");
			((KML_Unknown *)kmlobj->GetObjectByOffset(i))->DeepEcho(indent2);
            if(flag) _printf0_(indent << "        kmlobj: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "        kmlobj: [empty]\n");

	return;
}
/*}}}*/
void  KML_Object::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent2[81];

//  attributes always written in keyword line of derived classes
//  comments always written after keyword line of derived classes

/*  loop over the unknown objects for the object  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (kmlobj->Size())
		for (i=0; i<kmlobj->Size(); i++) {
			((KML_Unknown *)kmlobj->GetObjectByOffset(i))->Write(filout,indent2);
		}

	return;
}
/*}}}*/
void  KML_Object::Read(FILE* fid,char* kstr){/*{{{*/

	KML_Object*  kobj;

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</Object", 8))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_Object::Read -- Unexpected closing tag " << kstr << ".\n");}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_Object::Read -- Unexpected field \"" << kstr << "\".\n");}

	else if (!strncmp(kstr,"<Placemark",10)) {
		kobj=(KML_Object*)new KML_Placemark();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Folder", 7)) {
		kobj=(KML_Object*)new KML_Folder();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Document", 9)) {
		kobj=(KML_Object*)new KML_Document();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<GroundOverlay",14)) {
		kobj=(KML_Object*)new KML_GroundOverlay();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<LatLonBox",10)) {
		kobj=(KML_Object*)new KML_LatLonBox();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Icon", 5)) {
		kobj=(KML_Object*)new KML_Icon();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Point", 6)) {
		kobj=(KML_Object*)new KML_Point();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<LineString",11)) {
		kobj=(KML_Object*)new KML_LineString();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<LinearRing",11)) {
		kobj=(KML_Object*)new KML_LinearRing();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Polygon", 8)) {
		kobj=(KML_Object*)new KML_Polygon();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<MultiGeometry",14)) {
		kobj=(KML_Object*)new KML_MultiGeometry();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

//	else if (!strncmp(kstr,"<IconStyle",10)) {
//		kobj=(KML_Object*)new KML_IconStyle();
//		kobj->Read(fid,kstr);
//		kmlobj    ->AddObject((Object*)kobj);
//	}

//	else if (!strncmp(kstr,"<LabelStyle",11)) {
//		kobj=(KML_Object*)new KML_LabelStyle();
//		kobj->Read(fid,kstr);
//		kmlobj    ->AddObject((Object*)kobj);
//	}

	else if (!strncmp(kstr,"<LineStyle",10)) {
		kobj=(KML_Object*)new KML_LineStyle();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<PolyStyle",10)) {
		kobj=(KML_Object*)new KML_PolyStyle();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

//	else if (!strncmp(kstr,"<BalloonStyle",13)) {
//		kobj=(KML_Object*)new KML_BalloonStyle();
//		kobj->Read(fid,kstr);
//		kmlobj    ->AddObject((Object*)kobj);
//	}

//	else if (!strncmp(kstr,"<ListStyle",10)) {
//		kobj=(KML_Object*)new KML_ListStyle();
//		kobj->Read(fid,kstr);
//		kmlobj    ->AddObject((Object*)kobj);
//	}

	else if (!strncmp(kstr,"<",1)) {
		_printf0_("KML_Object::Read -- Unrecognized opening tag " << kstr << ".\n");
//		KMLFileTagSkip(kstr,
//					   fid);
		kobj=(KML_Object*)new KML_Unknown();
		kobj->Read(fid,kstr);
		kmlobj    ->AddObject((Object*)kobj);
	}

	return;
}
/*}}}*/
void  KML_Object::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	;

	return;
}
/*}}}*/
void  KML_Object::AddAttrib(const char* name,const char* value){/*{{{*/

	KML_Attribute* katt=NULL;

	katt=new KML_Attribute();
	katt->Alloc(name,value);
	katt->Add(attrib);

	return;
}
/*}}}*/
void  KML_Object::WriteAttrib(FILE* filout,const char* indent){/*{{{*/

//  attributes always written in keyword line of kml_object

/*  loop over any attributes for the object  */

	if (attrib->Size())
		for (int i=0; i<attrib->Size(); i++)
			((KML_Attribute *)attrib->GetObjectByOffset(i))->Write(filout,indent);

	return;
}
/*}}}*/
void  KML_Object::AddCommnt(int ncom,char** pcom){/*{{{*/

	int   i;
	KML_Comment* kcom=NULL;

	for (i=0; i<ncom; i++) {
		kcom=new KML_Comment();
		kcom->Alloc(pcom[i]);
		kcom->Add(commnt);
	}

	return;
}
/*}}}*/
void  KML_Object::AddCommnt(char* value){/*{{{*/

	KML_Comment* kcom=NULL;

	kcom=new KML_Comment();
	kcom->Alloc(value);
	kcom->Add(commnt);

	return;
}
/*}}}*/
void  KML_Object::WriteCommnt(FILE* filout,const char* indent){/*{{{*/

	int   i;

//  comments always written after keyword line of kml_object

/*  loop over any comments for the object  */

	if (commnt->Size())
		for (i=0; i<commnt->Size(); i++)
			((KML_Comment *)commnt->GetObjectByOffset(i))->Write(filout,indent);

	return;
}
/*}}}*/
