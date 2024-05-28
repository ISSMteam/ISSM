/*!\file KML_Comment.cpp
 * \brief: implementation of the kml_comment object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Comment.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Comment::KML_Comment(){/*{{{*/

	value     =NULL;

}
/*}}}*/
KML_Comment::~KML_Comment(){/*{{{*/

	if (value     ) xDelete<char>(value);

}
/*}}}*/

/*Other*/
void  KML_Comment::Echo(){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_("    ");
	if(flag) _printf0_(value << "\n");

	return;
}
/*}}}*/
void  KML_Comment::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Comment::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Comment::DeepEcho(const char* indent){/*{{{*/

	bool  flag=true;

	if(flag) _printf0_(indent << "    ");
	if(flag) _printf0_(value << "\n");

	return;
}
/*}}}*/
void  KML_Comment::Write(FILE* filout,const char* indent){/*{{{*/

	if (strncmp(&value[0]              ,"<!--",4))
		fprintf(filout,"%s<!--\n",indent);
	fprintf(filout,"%s  %s\n",indent,value);
	if (strncmp(&value[strlen(value)-3],"-->" ,3))
		fprintf(filout,"%s-->\n",indent);

	return;
}
/*}}}*/
void  KML_Comment::Read(FILE* fid,char* kstr){/*{{{*/

//  comments always read as part of KMLFileToken

	;

	return;
}
/*}}}*/
void  KML_Comment::Alloc(const char* valuei){/*{{{*/

	value=xNew<char>(strlen(valuei)+1);
	memcpy(value,valuei,(strlen(valuei)+1)*sizeof(char));

	return;
}
/*}}}*/
void  KML_Comment::Add(DataSet* commnt){/*{{{*/

	commnt->AddObject((Object*)this);

	return;
}
/*}}}*/
void  KML_Comment::Get(char** pvalueo){/*{{{*/

	*pvalueo=xNew<char>(strlen(value)+1);
	memcpy(*pvalueo,value,(strlen(value)+1)*sizeof(char));

	return;
}
/*}}}*/
