/*!\file KML_Unknown.cpp
 * \brief: implementation of the kml_unknown object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KMLFileReadUtils.h"
#include "./KML_Unknown.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Unknown::KML_Unknown(){/*{{{*/

	name      =NULL;
	value     =NULL;

}
/*}}}*/
KML_Unknown::~KML_Unknown(){/*{{{*/

	if (name      ) xDelete<char>(name);
	if (value     ) xDelete<char>(value);

}
/*}}}*/

/*Other*/
void  KML_Unknown::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("KML_Unknown " << name << ":\n");
	KML_Object::Echo();

	if(value){
		if(flag) _printf0_("         value: \"" << value << "\"\n");
	}
	else{
		if(flag) _printf0_("         value: [none]\n");
	}

	return;
}
/*}}}*/
void  KML_Unknown::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Unknown::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Unknown::DeepEcho(const char* indent){/*{{{*/

	char*        valuei;
	char*        vtoken;
	char         nl[]={'\n','\0'};
	bool         flag=true;

	if(flag) _printf0_(indent << "KML_Unknown " << name << ":\n");
	KML_Object::DeepEcho(indent);

	if (value     ) {
		valuei=xNew<char>(strlen(value)+1);
		memcpy(valuei,value,(strlen(value)+1)*sizeof(char)); 

		vtoken=strtok(valuei,nl);
		if(flag) _printf0_(indent << "         value: \"" << vtoken);

		while((vtoken=strtok(NULL,nl)))
			if(flag) _printf0_("\n" << indent << "                 " << vtoken);
		if(flag) _printf0_("\"\n");

		xDelete<char>(valuei);
	}
    else
        if(flag) _printf0_(indent << "         value: [none]\n");

	return;
}
/*}}}*/
void  KML_Unknown::Write(FILE* filout,const char* indent){/*{{{*/

	char*        valuei;
	char*        vtoken;
	char         nl[]={'\n','\0'};

	fprintf(filout,"%s<%s",indent,name);
	WriteAttrib(filout," ");
	fprintf(filout,">\n");
	WriteCommnt(filout,indent);

	if (value     ) {
		valuei=xNew<char>(strlen(value)+1);
		memcpy(valuei,value,(strlen(value)+1)*sizeof(char)); 

		vtoken=strtok(valuei,nl);
		fprintf(filout,"%s  %s\n",indent,vtoken);

		while((vtoken=strtok(NULL,nl)))
			fprintf(filout,"%s  %s\n",indent,vtoken);

		xDelete<char>(valuei);
	}

	KML_Object::Write(filout,indent);

	fprintf(filout,"%s</%s>\n",indent,name);

	return;
}
/*}}}*/
void  KML_Unknown::Read(FILE* fid,char* kstr){/*{{{*/

	char*        kstri;
	char*        value2=NULL;
	int          ncom=0;
	char**       pcom=NULL;
	char         nl[]={'\n','\0'};

/*  get object name  */

	name=KMLFileTagName(NULL,
						kstr);
//	_printf0_("KML_Unknown::Read -- opening name=" << name << ".\n");

/*  get object attributes and check for solo tag  */

	if (KMLFileTagAttrib(this,
						 kstr))
		return;

/*  loop over and process fields within opening and closing tags  */

	while((kstri=KMLFileToken(fid, &ncom,&pcom))){
//		_printf0_("KML_Unknown::Read -- kstri=" << kstri << ".\n");
		if      (!strncmp(&kstri[0],"</", 2) &&
				 !strncmp(&kstri[2],name,strlen(name))) {
//			_printf0_("KML_Unknown::Read -- closing name=" << name << ".\n");
			xDelete<char>(kstri);
			break;
		}
		else if (!strncmp(kstri,"</",2))
		  {_error_("KML_Unknown::Read -- Unexpected closing tag " << kstri << ".\n");}

		else if (strncmp(kstri,"<",1)) {
			if (value) {
				value2=xNew<char>(strlen(value)+1+strlen(kstri)+1);
				memcpy(value2,value,(strlen(value)+1)*sizeof(char));
				xDelete<char>(value);
				value=value2;
				value2=NULL;
//				value=(char *) xrealloc(value,(strlen(value)+1+strlen(kstri)+1)*sizeof(char));
				strcat(value,nl);
				strcat(value,kstri);
			}
			else {
				value=xNew<char>(strlen(kstri)+1);
				memcpy(value,kstri,(strlen(kstri)+1)*sizeof(char));
			}
		}

		else if (!strncmp(kstri,"<",1))
			KML_Object::Read(fid,kstri);

		xDelete<char>(kstri);
	}

	this->AddCommnt(ncom,pcom);

	for(ncom=ncom; ncom>0; ncom--)
		xDelete<char>(pcom[ncom-1]);
	xDelete<char*>(pcom);

	return;
}
/*}}}*/
