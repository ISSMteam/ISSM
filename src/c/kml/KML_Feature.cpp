/*!\file KML_Feature.cpp
 * \brief: implementation of the kml_feature abstract object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../shared/shared.h"
/*}}}*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./KML_Feature.h"
#include "./KML_Style.h"
#include "./KMLFileReadUtils.h"
#include "../shared/shared.h"
/*}}}*/

/*Constructors/destructor/copy*/
KML_Feature::KML_Feature(){/*{{{*/

	memcpy(name,"",(strlen("")+1)*sizeof(char));

	visibility=true;
	open      =false;
	memcpy(snippet,"",(strlen("")+1)*sizeof(char));
	memcpy(descript,"",(strlen("")+1)*sizeof(char));
	memcpy(styleurl,"",(strlen("")+1)*sizeof(char));
	style     =new DataSet;

}
/*}}}*/
KML_Feature::~KML_Feature(){/*{{{*/

	if (style) {
		delete style;
		style     =NULL;
	}

}
/*}}}*/

/*Other*/
void  KML_Feature::Echo(){/*{{{*/

	bool  flag=true;

	KML_Object::Echo();

	if(flag) _printf0_("          name: \"" << name << "\"\n");
	if(flag) _printf0_("    visibility: " << (visibility ? "true" : "false") << "\n");
	if(flag) _printf0_("          open: " << (open ? "true" : "false") << "\n");
	if(flag) _printf0_("       snippet: \"" << snippet << "\"\n");
	if(flag) _printf0_("      descript: \"" << descript << "\"\n");
	if(flag) _printf0_("      styleurl: \"" << styleurl << "\"\n");
	if(flag) _printf0_("         style: (size=" << style->Size() << ")\n");

	return;
}
/*}}}*/
void  KML_Feature::DeepEcho(){/*{{{*/

	char  indent[81]="";

	KML_Feature::DeepEcho(indent);

	return;
}
/*}}}*/
void  KML_Feature::DeepEcho(const char* indent){/*{{{*/

	int   i;
	char  indent2[81];
	bool  flag=true;

	KML_Object::DeepEcho(indent);

	if(flag) _printf0_(indent << "          name: \"" << name << "\"\n");
	if(flag) _printf0_(indent << "    visibility: " << (visibility ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "          open: " << (open ? "true" : "false") << "\n");
	if(flag) _printf0_(indent << "       snippet: \"" << snippet << "\"\n");
	if(flag) _printf0_(indent << "      descript: \"" << descript << "\"\n");
	if(flag) _printf0_(indent << "      styleurl: \"" << styleurl << "\"\n");

/*  loop over any styles for the feature  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));
	strcat(indent2,"  ");

	if (style->Size())
		for (i=0; i<style->Size(); i++) {
			if(flag) _printf0_(indent << "         style: -------- begin [" << i << "] --------\n");
			((KML_Style *)style->GetObjectByOffset(i))->DeepEcho(indent2);
			if(flag) _printf0_(indent << "         style: --------  end  [" << i << "] --------\n");
		}
	else
		if(flag) _printf0_(indent << "         style: [empty]\n");

	return;
}
/*}}}*/
void  KML_Feature::Write(FILE* filout,const char* indent){/*{{{*/

	int   i;
	char  indent2[81];

	KML_Object::Write(filout,indent);

	if (name     && strlen(name))
		fprintf(filout,"%s  <name>%s</name>\n",indent,name);
	fprintf(filout,"%s  <visibility>%d</visibility>\n",indent,(visibility ? 1 : 0));
	fprintf(filout,"%s  <open>%d</open>\n",indent,(open ? 1 : 0));
	if (snippet  && strlen(snippet))
		fprintf(filout,"%s  <Snippet maxLines=\"2\">%s</Snippet>\n",indent,snippet);
	if (descript && strlen(descript))
		fprintf(filout,"%s  <description>%s</description>\n",indent,descript);
	if (styleurl && strlen(styleurl))
		fprintf(filout,"%s  <styleUrl>%s</styleUrl>\n",indent,styleurl);

/*  loop over any styles for the feature  */

	memcpy(indent2,indent,(strlen(indent)+1)*sizeof(char));

	strcat(indent2,"  ");

    for (i=0; i<style->Size(); i++)
        ((KML_Style *)style->GetObjectByOffset(i))->Write(filout,indent2);

	return;
}
/*}}}*/
void  KML_Feature::Read(FILE* fid,char* kstr){/*{{{*/

	KML_Object*  kobj;

/*  process field within opening and closing tags  */

	if      (!strncmp(kstr,"</Feature", 9))
		return;
	else if (!strncmp(kstr,"</",2))
	  {_error_("KML_Feature::Read -- Unexpected closing tag " << kstr);}
	else if (strncmp(kstr,"<",1))
	  {_error_("KML_Feature::Read -- Unexpected field \"" << kstr << "\"");}

	else if (!strncmp(kstr,"<Style", 6)) {
		kobj=(KML_Object*)new KML_Style();
		kobj->Read(fid,kstr);
		style     ->AddObject((Object*)kobj);
	}

	else if (!strcmp(kstr,"<name>"))
		KMLFileTokenParse( name      ,NULL,KML_FEATURE_NAME_LENGTH, kstr, fid);
	else if (!strcmp(kstr,"<visibility>"))
		KMLFileTokenParse(&visibility, kstr, fid);
	else if (!strcmp(kstr,"<open>"))
		KMLFileTokenParse(&open      , kstr, fid);
	else if (!strncmp(kstr,"<snippet", 8))
		KMLFileTokenParse( snippet   ,NULL,KML_FEATURE_SNIPPET_LENGTH, kstr, fid);
	else if (!strcmp(kstr,"<description>"))
		KMLFileTokenParse( descript  ,NULL,KML_FEATURE_DESCRIPT_LENGTH, kstr, fid);
	else if (!strcmp(kstr,"<styleUrl>"))
		KMLFileTokenParse( styleurl  ,NULL,KML_FEATURE_STYLEURL_LENGTH, kstr, fid);

	else if (!strncmp(kstr,"<",1))
		KML_Object::Read(fid,kstr);

	return;
}
/*}}}*/
