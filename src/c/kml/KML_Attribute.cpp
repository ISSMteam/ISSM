/*!\file KML_Attribute.cpp
 * \brief: implementation of the kml_attribute object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Attribute.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Attribute::KML_Attribute(){/*{{{*/

	name      =NULL;
	value     =NULL;

}
/*}}}*/
KML_Attribute::~KML_Attribute(){/*{{{*/

	if (name      ) xDelete<char>(name);
	if (value     ) xDelete<char>(value);

}
/*}}}*/

/*Other*/
void  KML_Attribute::Echo(){/*{{{*/

	int   i;
	bool  flag=true;

	if(flag) _printf0_("    ");
	for (i=0;i<10-strlen(name);i++)
		if(flag) _printf0_(" ");
	if(flag) _printf0_(name << ": \"" << value << "\"\n");

	return;
}
/*}}}*/
void  KML_Attribute::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Attribute::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Attribute::DeepEcho(const char* indent){/*{{{*/

	int   i;
	bool  flag=true;

	if(flag) _printf0_(indent << "    ");
	for (i=0;i<10-strlen(name);i++)
		if(flag) _printf0_(" ");
	if(flag) _printf0_(name << ": \"" << value << "\"\n");

	return;
}
/*}}}*/
void  KML_Attribute::Write(FILE* filout,const char* indent){/*{{{*/

//  attributes always written in keyword line of kml_object

	fprintf(filout,"%s%s=\"%s\"",indent,name,value);

	return;
}
/*}}}*/
void  KML_Attribute::Read(FILE* fid,char* kstr){/*{{{*/

//  attributes always read in keyword line of kml_object

	;

	return;
}
/*}}}*/
void  KML_Attribute::Alloc(const char* namei,const char* valuei){/*{{{*/

	name =xNew<char>(strlen(namei )+1);
	memcpy(name,namei,(strlen(namei)+1)*sizeof(char));

	value=xNew<char>(strlen(valuei)+1);
	memcpy(value,valuei,(strlen(valuei)+1)*sizeof(char));

	return;
}
/*}}}*/
void  KML_Attribute::Add(DataSet* attrib){/*{{{*/

	attrib->AddObject((Object*)this);

	return;
}
/*}}}*/
void  KML_Attribute::Get(char** pvalueo,char* deflt){/*{{{*/

	if (!value || !strlen(value)) {
		*pvalueo=xNew<char>(strlen(deflt)+1);
		memcpy(*pvalueo,deflt,(strlen(deflt)+1)*sizeof(char));
	}
	else {
		*pvalueo=xNew<char>(strlen(value)+1);
		memcpy(*pvalueo,value,(strlen(value)+1)*sizeof(char));
	}

	return;
}
/*}}}*/
