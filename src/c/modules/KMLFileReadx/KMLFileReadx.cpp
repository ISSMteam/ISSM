/*!\file KMLFileReadx.cpp
 */

#include "./KMLFileReadx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

KML_Object* KMLFileReadx(FILE* fid){

	char*   kstr;
	KML_File*      kxml=NULL;
	KML_File*      kdtd=NULL;
	KML_File*      kfil=NULL;

	clock_t clock0,clock1;
	time_t  time0, time1;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nKMLFileReadx Module -- " << ctime(&time0));

/*  read kml file  */

	while((kstr=KMLFileToken(fid, NULL,NULL))){
		if      (!strncmp(kstr,"<?xml"    ,5)) {
			kxml=new KML_File();
			KMLFileTagAttrib(kxml,
							 kstr);
		}
		else if (!strncmp(kstr,"<!DOCTYPE",9)) {
			kdtd=new KML_File();
			KMLFileTagAttrib(kdtd,
							 kstr);
		}
		else if (!strncmp(kstr,"<kml"     ,4)) {
			kfil=new KML_File();
			kfil->Read(fid,kstr);
//			kfil->DeepEcho();
		}

//		_printf0_(kstr << "\n");
		xDelete<char>(kstr);
	}

	if (kxml) {
		_printf0_("XML declaration:\n");
		kxml->DeepEcho("  ");
		delete kxml;
	}
	if (kdtd) {
		_printf0_("DTD declaration (not yet implemented):\n");
		kdtd->DeepEcho("  ");
		delete kdtd;
	}

	clock1=clock();
	time1 =time(NULL);
	_printf_("KMLFileReadx Module -- " <<((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " <<difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return(kfil);
}
