/*!\file KML_Container.cpp
 * \brief: implementation of the kml_container abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Placemark.h"
#include "./KML_Document.h"
#include "./KML_GroundOverlay.h"
#include "./KML_Folder.h"
#include "./KML_Object.h"
#include "./KML_Container.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Container::KML_Container(){/*{{{*/

	feature   =new DataSet;

}
/*}}}*/
KML_Container::~KML_Container(){/*{{{*/

	if (feature) {
		delete feature;
		feature   =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Container::Echo(){/*{{{*/

	bool  flag=true;

	KML_Feature::Echo();

	if(flag) _printf0_("       feature: (size=" << feature->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_Container::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Container::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Container::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

	KML_Feature::DeepEcho(indent);

/*  loop over the features for the container  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (feature->Size())
		for (i=0; i<feature->Size(); i++) {
			if(flag) _printf0_(indent << "       feature: -------- begin [" << i << "] --------\n");
			((KML_Feature *)feature->GetObjectByOffset(i))->DeepEcho(indent2);
			if(flag) _printf0_(indent << "       feature: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "       feature: [empty]\n");

	return;
}
/*}}}*/
void  KML_Container::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent2[81];

	KML_Feature::Write(filout,indent);

/*  loop over the features for the container  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

	for (i=0; i<feature->Size(); i++)
		((KML_Feature *)feature->GetObjectByOffset(i))->Write(filout,indent2);

	return;
}
/*}}}*/
void  KML_Container::Read(FILE* fid,char* kstr){/*{{{*/

	KML_Object*  kobj;

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</Container",11)) {
		xDelete<char>(kstr);
		return;
	}
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_Container::Read -- Unexpected closing tag " << kstr );}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_Container::Read -- Unexpected field \"" << kstr << "\"");}

	else if (!strncmp(kstr,"<Placemark",10)) {
		kobj=(KML_Object*)new KML_Placemark();
		kobj->Read(fid,kstr);
		feature   ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Folder", 7)) {
		kobj=(KML_Object*)new KML_Folder();
		kobj->Read(fid,kstr);
		feature   ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<Document", 9)) {
		kobj=(KML_Object*)new KML_Document();
		kobj->Read(fid,kstr);
		feature   ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<GroundOverlay",14)) {
		kobj=(KML_Object*)new KML_GroundOverlay();
		kobj->Read(fid,kstr);
		feature   ->AddObject((Object*)kobj);
	}

	else if (!strncmp(kstr,"<",1))
		KML_Feature::Read(fid,kstr);

	return;
}
/*}}}*/
void  KML_Container::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int   i;

/*  loop over the features for the container  */

	for (i=0; i<feature->Size(); i++)
		((KML_Object *)feature->GetObjectByOffset(i))->WriteExp(fid,nstr,sgn,cm,sp);

	return;
}
/*}}}*/
