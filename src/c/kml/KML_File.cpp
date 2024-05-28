/*!\file KML_File.cpp
 * \brief: implementation of the kml_file object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_File.h" 
#include "./KML_Object.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_File::KML_File(){/*{{{*/

	;

}
/*}}}*/
KML_File::~KML_File(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_File::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_File:\n");
	KML_Object::Echo();

	return;
}
/*}}}*/
void  KML_File::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_File::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_File::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_File:\n");
	KML_Object::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_File::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<kml",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Object::Write(filout,indent);

	fprintf(filout,"%s</kml>\n",indent);

	return;
}
/*}}}*/
void  KML_File::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</kml", 5)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_File::Read -- Unexpected closing tag " << kstri << ".");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_File::Read -- Unexpected field \"" << kstri << "\"");}

		else if (!strncmp(kstri,"<",1))
			KML_Object::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(; ncom>0; ncom--) xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
void  KML_File::WriteExp(FILE* fid,const char* nstr,int sgn,double cm,double sp){/*{{{*/

	int   i;

/*  loop over the kml objects for the file  */

	for (i=0; i<kmlobj->Size(); i++)
		((KML_Object *)kmlobj->GetObjectByOffset(i))->WriteExp(fid,nstr,sgn,cm,sp);

	return;
}
/*}}}*/
