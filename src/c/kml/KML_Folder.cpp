/*!\file KML_Folder.cpp
 * \brief: implementation of the kml_folder object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Folder.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Folder::KML_Folder(){/*{{{*/

	;

}
/*}}}*/
KML_Folder::~KML_Folder(){/*{{{*/

	;

}
/*}}}*/

/*Other*/
void  KML_Folder::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Folder:\n");
	KML_Container::Echo();

	return;
}
/*}}}*/
void  KML_Folder::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Folder::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Folder::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "KML_Folder:\n");
	KML_Container::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Folder::Write(FILE* filout,const char* indent){/*{{{*/

	fprintf(filout,"%s<Folder",indent);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	KML_Container::Write(filout,indent);

	fprintf(filout,"%s</Folder>\n",indent);

	return;
}
/*}}}*/
void  KML_Folder::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	int          ncom=0;
	char**       pcom=NULL;

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
		if      (!strncmp(kstri,"</Folder", 8)) {
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Folder::Read -- Unexpected closing tag " << kstri << ".\n");}
		else if (strncmp(kstri,"<",1))
		  {_error_("KML_Folder::Read -- Unexpected field \"" << kstri << "\".\n");}

		else if (!strncmp(kstri,"<",1))
			KML_Container::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for (ncom=ncom; ncom>0; ncom--) xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
